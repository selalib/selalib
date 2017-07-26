!> @ingroup splines
!> Implements arbitrary degree bspline interpolation on a uniform grid
!> given a B-Spline object from sll_m_bsplines
module sll_m_spline_1d_non_uniform

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
  sll_p_periodic, &
  sll_p_hermite , &
  sll_p_greville, &
  sll_p_open    , &
  sll_p_mirror

use sll_m_spline_1d_base, only: &
  sll_c_spline_1d

use sll_m_bsplines, only: &
  sll_t_bsplines, &
  sll_s_bsplines_init_from_grid, &
  sll_s_bsplines_free, &
  sll_f_find_cell, &
  sll_s_bsplines_eval_basis, &
  sll_s_bsplines_eval_deriv, &
  sll_s_bsplines_eval_basis_and_n_derivs, &
  sll_s_uniform_bsplines_eval_basis

use schur_complement, only: &
  schur_complement_solver, &
  schur_complement_fac   , &
  schur_complement_slv   , &
  schur_complement_free

implicit none

public :: &
  sll_t_spline_1d_non_uniform

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Allowed boundary conditions
  sll_int32, parameter :: &
    allowed_bcs(*) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

!> @brief
!> basic type for one-dimensional B-spline interpolation.
type, extends(sll_c_spline_1d) :: sll_t_spline_1d_non_uniform

  type (sll_t_bsplines) :: bsp
  sll_int32                               :: n        ! dimension of spline space
  sll_int32                               :: deg      ! degree of spline
  sll_real64, allocatable                 :: tau(:)   ! interpolation points
  sll_real64, allocatable                 :: bcoef(:) ! bspline coefficients
  sll_int32                               :: bc_type  ! boundary condition
  sll_real64, allocatable                 :: q(:,:)   ! triangular factorization
  sll_real64, allocatable                 :: bsdx(:,:)
  sll_int32,  allocatable                 :: ipiv(:)
  sll_int32                               :: offset
  type(schur_complement_solver)           :: schur

contains

  procedure :: init                => s_spline_1d_non_uniform__init
  procedure :: free                => s_spline_1d_non_uniform__free
  procedure :: compute_interpolant => s_spline_1d_non_uniform__compute_interpolant
  procedure :: eval                => f_spline_1d_non_uniform__eval
  procedure :: eval_deriv          => f_spline_1d_non_uniform__eval_deriv
  procedure :: eval_array          => s_spline_1d_non_uniform__eval_array
  procedure :: eval_array_deriv    => s_spline_1d_non_uniform__eval_array_deriv
  procedure :: get_coeff           => f_spline_1d_non_uniform__get_coeff
  procedure :: get_interp_points   => s_spline_1d_non_uniform__get_interp_points

end type sll_t_spline_1d_non_uniform

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function f_spline_1d_non_uniform__get_coeff( self ) result( ptr )

    class(sll_t_spline_1d_non_uniform), target, intent(in) :: self
    sll_real64, pointer :: ptr(:)

    ptr => self%bcoef

  end function f_spline_1d_non_uniform__get_coeff

  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_non_uniform__get_interp_points( self, tau )

    class(sll_t_spline_1d_non_uniform), intent(in   ) :: self
    sll_real64,            allocatable, intent(  out) :: tau(:)

    SLL_ASSERT( allocated( self%tau ) )
    allocate( tau(size(self%tau)), source=self%tau )

  end subroutine s_spline_1d_non_uniform__get_interp_points

  !-----------------------------------------------------------------------------
  !> @brief     Constructor for sll_t_bspline_1d object
  !> @param[in] num_pts  Number of points where the data to be
  !>                     interpolated are represented.
  !> @param[in] degree   Spline degree
  !> @param[in] xmin     Minimum value of the abscissae where data are meant
  !>                     to be interpolated.
  !> @param[in] xmax     Maximum value of the abscissae where data are meant
  !>                     to be interpolated.
  !> @param[in] bc_xmin  Boundary condition at x=xmin
  !> @param[in] bc_xmax  Boundary condition at x=xmax
  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_non_uniform__init( &
    self   , &
    degree , &
    breaks , &
    bc_xmin, &
    bc_xmax )

    class(sll_t_spline_1d_non_uniform), intent(  out) :: self
    sll_int32                         , intent(in   ) :: degree
    sll_real64                        , intent(in   ) :: breaks(:)
    sll_int32                         , intent(in   ) :: bc_xmin
    sll_int32                         , intent(in   ) :: bc_xmax

    sll_int32  :: num_pts
    sll_int32  :: i, j
!    sll_real64 :: delta
    sll_int32  :: bc_type
    sll_int32  :: basis_bc_xmin
    sll_int32  :: basis_bc_xmax

    ! Sanity checks
    SLL_ASSERT( degree >= 1 )
    SLL_ASSERT( size( breaks ) >= 2 )
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )

    ! NOTE: in the future different boundary conditions at xmin and xmax
    !       should be considered. For now we only check that bc_xmin==bc_xmax
    SLL_ASSERT( bc_xmin == bc_xmax )
    bc_type = bc_xmin

    num_pts = size( breaks )

    ! set first attributes
    self%bc_type = bc_type
    if (bc_type == sll_p_periodic) then
       self%n = num_pts - 1 ! dimension of periodic spline space
       self%offset = degree/2  ! offset needed for periodic spline evaluation
    else
       self%n = num_pts + degree - 1 ! dimension of non periodic spline space
       self%offset = 0
    end if
    self%deg = degree

    ! Determine boundary conditions for B-splines
    select case (bc_xmin)
    case (sll_p_periodic); basis_bc_xmin = sll_p_periodic
    case (sll_p_greville); basis_bc_xmin = sll_p_open
    case (sll_p_hermite ); basis_bc_xmin = sll_p_open ! sll_p_mirror also possible
    end select

    select case (bc_xmax)
    case (sll_p_periodic); basis_bc_xmax = sll_p_periodic
    case (sll_p_greville); basis_bc_xmax = sll_p_open
    case (sll_p_hermite ); basis_bc_xmax = sll_p_open ! sll_p_mirror also possible
    end select

    ! Construct a sll_t_bsplines object
    call sll_s_bsplines_init_from_grid( self%bsp, degree, breaks,  &
          basis_bc_xmin, basis_bc_xmax )

    allocate( self%bcoef(self%n) )
    allocate( self%bsdx (self%deg/2+1,self%deg+1) )
    allocate( self%ipiv (self%n) )

    ! Determine array 'tau' of interpolation points
    select case (bc_type)

    case (sll_p_periodic)

      ! Even degree: interpolation points are cell midpoints
      ! Odd  degree: interpolation points are grid points excluding last
      ! TODO: use averaging (Greville style)
      allocate( self%tau(self%n) )
      if (modulo(degree,2) == 0) then
        self%tau(:) = [(0.5_f64*(breaks(i)+breaks(i+1)), i=1, num_pts-1)]
      else
        self%tau(:) = breaks(1:num_pts-1)
      end if

    case (sll_p_hermite)

       ! In this case only interior points are saved in tau
       ! Even degree: interpolation points are cell midpoints
       ! Odd  degree: interpolation points are grid points including last
       if (modulo(degree,2) == 0) then
         allocate( self%tau(num_pts-1) )
         self%tau(:) = [(0.5_f64*(breaks(i)+breaks(i+1)), i=1, num_pts-1)]
       else
         allocate( self%tau(num_pts) )
         self%tau(:) = breaks(1:num_pts)
       end if

    case (sll_p_greville)

      ! Set interpolation points to be Greville points
      allocate( self%tau(self%n) )
      do i = 1, self%n
        self%tau(i) = sum( self%bsp%knots(i+1-degree:i) ) / real(degree,f64)
      end do

    end select

    ! Special case: linear spline
    ! No need for matrix assembly
    if (self%deg == 1) then
      allocate( self%q(0,0) )
      return
    end if

    ! Assemble banded matrix (B_j(tau(i))) for spline interpolation
    select case (bc_type)
    case (sll_p_periodic)
       call build_system_periodic(self)
    case (sll_p_greville)
       call build_system_greville(self)
    case (sll_p_hermite)
       call build_system_with_derivative(self)
    end select

  end subroutine s_spline_1d_non_uniform__init

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system needed for periodic spline interpolation
  !> @param[inout] self bspline interpolation object
  !-----------------------------------------------------------------------------
  subroutine build_system_periodic( self )

    class(sll_t_spline_1d_non_uniform), intent(inout) :: self

    integer  :: i, j
    integer  :: k, s
    integer  :: icell
    real(wp) :: x
    real(wp) :: values(self%deg+1)

    k = self%deg/2
    allocate( self%q(2*k+1,self%n) )

    ! evaluate bsplines at interpolation points
    ! and assemble matrix q for linear system
    do i = 1, self%n
      x = self%tau(i)
      icell = sll_f_find_cell( self%bsp, x )
      call sll_s_bsplines_eval_basis(self%bsp, icell, x, values)
      do s = 1, 2*k+1
        j = modulo(icell-k-2+s,self%n)+1
        self%q(2*k+2-s,j) = values(s)
      end do
    end do

    ! Perform LU decomposition of matrix q
    call schur_complement_fac( self%schur, self%n, k, self%q )

  end subroutine build_system_periodic

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for spline interpolation at Greville points
  !> @param[inout] self bspline object
  !-----------------------------------------------------------------------------
  subroutine build_system_greville( self )

    class(sll_t_spline_1d_non_uniform), intent(inout) :: self

    sll_int32  :: iflag
    sll_int32  :: j
    sll_int32  :: k
    sll_int32  :: ii
    sll_int32  :: jj
    sll_int32  :: icell
    sll_real64 :: x
    sll_real64 :: values(self%deg+1)

    ! The matrix q is a banded matrix using the storage required by banfac (De Boor)
    ! It has k bands above diagonal, k bands below and the diagonal itself
    ! The term A(ii,jj) of the full matrix is stored in q(ii-jj+k+1,jj)
    ! The Bspline interpolation matrix at Greville points x_ii is
    ! A(ii,jj) = B_jj(x_ii)
    k = self%deg - 1
    allocate( self%q(2*k+1,self%n) )
    self%q = 0.0_f64

    ! Treat i=1 separately
    x =  self%bsp%xmin
    icell = 1
    ii = 1   ! first Greville point
    call sll_s_bsplines_eval_basis(self%bsp, icell, x, values)
    ! iterate only to k+1 as last bspline is 0
    do j=1,k+1
       jj = icell+j-1
       self%q(ii-jj+k+1,jj) = values(j)
    end do

    do ii=2,self%n - 1
       x = self%tau(ii)
       icell = sll_f_find_cell(self%bsp, x )
       call sll_s_bsplines_eval_basis(self%bsp, icell, x, values)
       do j=1,k+2
          jj = icell+j-1
          self%q(ii-jj+k+1,jj) = values(j)
       end do
    end do
    ! Treat i=self%n separately
    x =  self%bsp%xmax
    icell = self%n - self%deg
    ii = self%n  ! last Greville point
    call sll_s_bsplines_eval_basis(self%bsp, icell, x, values)
    ! iterate only to k as first bspline is 0
    do j=2,k+2
       jj = icell+j-1
       self%q(ii-jj+k+1,jj) = values(j)
    end do
    ! Perform LU decomposition of matrix q
    call banfac ( self%q, 2*k+1, self%n, k, k, iflag )

  end subroutine build_system_greville

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for spline interpolation with Hermite boundary conditions
  !> @param[inout] self bspline object
  !-----------------------------------------------------------------------------
  subroutine build_system_with_derivative( self )

    class(sll_t_spline_1d_non_uniform), intent(inout) :: self

    sll_int32  :: nbc
    sll_int32  :: iflag
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: ii
    sll_int32  :: jj
    sll_int32  :: offset
    sll_int32  :: k
    sll_int32  :: icell
    sll_real64 :: x
    sll_real64 :: values(self%deg+1)

    ! number of boundary conditions needed depending on spline degree
    k   = self%deg-1
    nbc = self%deg/2

    ! The matrix q is a banded matrix using the storage required by DGBTRF (LAPACK)
    ! It has 2*k bands above diagonal, k bands below and the diagonal itself
    ! k additional bands above diagonal needed for pivoting
    ! The term A(ii,jj) of the full matrix is stored in q(ii-jj+2*k+1,jj)
    ! The Bspline interpolation matrix at Greville points x_ii is
    ! A(ii,jj) = B_jj(x_ii)
    allocate( self%q(3*k+1,self%n) )
    self%q = 0.0_f64

    ! For even degree splines interpolation points are at cell midpoints
    ! value of function on the boundary needed as additional boundary conditions
    ! for odd degree splines  baundary is in the interpolation points
    ! only derivative values are needed as boundary conditions
    if (modulo(k+1,2)==0) then ! spline degree even
       offset = 0
    else
       offset = 1
    end if

    ! boundary conditions at xmin
    ii=0
    x = self%bsp%xmin
    icell = 1
    call sll_s_bsplines_eval_basis_and_n_derivs( self%bsp, icell, x , nbc , self%bsdx )
    do i=1, nbc
       ! iterate only to k+1 as last bspline is 0
       ii=ii+1
       do j=1,k+1
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = self%bsdx(i+offset,j)
       end do
    end do
    ! interpolation points
    do i=1,self%n - 2*nbc
       ii = ii + 1
       x = self%tau(i)
       icell = sll_f_find_cell( self%bsp, x )
       call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
       do j=1,k+2
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = values(j)
       end do
    end do
    ! boundary conditions at xmax
    x = self%bsp%xmax
    icell = self%n - self%deg
    call sll_s_bsplines_eval_basis_and_n_derivs( self%bsp, icell, x , nbc , self%bsdx )
    do i= 1, nbc
       ii = ii + 1
       do j=2,k+2
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = self%bsdx(i+offset,j)
       end do
    end do

    ! Perform LU decomposition of matrix q with Lapack
    call dgbtrf( self%n, self%n, k, k, self%q, 3*k+1, self%ipiv, iflag )

  end subroutine build_system_with_derivative

  !-----------------------------------------------------------------------------
  !> @brief        Compute interpolating 1D bspline
  !> Computes coefficients of 1D bspline that interpolates function on grid.
  !> If Hermite BCs are used, derivatives at boundary are also needed.
  !>
  !> @param[inout] self       1D bspline object
  !> @param[in]    gtau       function values of interpolation points
  !> @param[in]    derivs_xmin (optional) array with boundary conditions at xmin
  !> @param[in]    derivs_xmax (optional) array with boundary conditions at xmax
  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_non_uniform__compute_interpolant( self, gtau, derivs_xmin, derivs_xmax )

    class(sll_t_spline_1d_non_uniform), intent(inout) :: self
    sll_real64                        , intent(in   ) :: gtau(:)
    sll_real64,               optional, intent(in   ) :: derivs_xmin(:)
    sll_real64,               optional, intent(in   ) :: derivs_xmax(:)

    sll_int32 :: ncond
    sll_int32 :: k
    sll_int32 :: iflag

    ! Special case: linear spline
    if (self%deg == 1) then
      self%bcoef(:) = gtau(1:self%n)
      return
    end if

    ! Degree > 1: boundary conditions matter
    select case (self%bc_type)
    case(sll_p_periodic)
       self%bcoef = gtau(1:self%n)
       call schur_complement_slv( self%schur, self%n, self%deg/2, self%q,  self%bcoef )
    case (sll_p_greville)
       self%bcoef = gtau
       call banslv( self%q, 2*self%deg-1, self%n, self%deg-1, self%deg-1, self%bcoef )
    case (sll_p_hermite)
       ! number of needed conditions at boundary
       ncond = self%deg/2
       if (present(derivs_xmin)) then
          self%bcoef(1:ncond) = derivs_xmin(1:ncond)
       else  ! set needed boundary values to 0
          self%bcoef(1:ncond) = 0.0_f64
       end if
       self%bcoef(ncond+1:self%n-ncond) = gtau(1:self%n-2*ncond)
       if (present(derivs_xmax)) then
          self%bcoef(self%n-ncond+1:self%n) = derivs_xmax(1:ncond)
       else ! set needed boundary values to 0
          self%bcoef(self%n-ncond+1:self%n) = 0.0_f64
       end if
       ! Use Lapack to solve banded system
       k = self%deg-1
       call dgbtrs( 'N', self%n, k, k, 1, self%q, 3*k+1, self%ipiv, &
                    self%bcoef, self%n, iflag )
    end select

  end subroutine s_spline_1d_non_uniform__compute_interpolant

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate spline S at given point x
  !> @param[in]  self spline object
  !> @param[in]  x    evaluation point
  !> @return     y    spline value y=S(x)
  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d_non_uniform__eval( self, x ) result( y )

    class(sll_t_spline_1d_non_uniform), intent(in) :: self
    sll_real64            , intent(in) :: x
    sll_real64 :: y

    sll_int32  :: icell
    sll_int32  :: ib
    sll_int32  :: j
    sll_real64 :: values(self%deg+1)

    icell = sll_f_find_cell( self%bsp, x )
    call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
    y = 0.0_f64
    do j = 1, self%deg+1
       ib = mod( icell+j-2-self%offset+self%n, self%n ) + 1
       y = y + values(j)*self%bcoef(ib)
    end do

  end function f_spline_1d_non_uniform__eval

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate spline derivative S' at given point x
  !> @param[in]  self spline object
  !> @param[in]  x    evaluation point
  !> @param[out] y    spline derivative y=S'(x)
  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d_non_uniform__eval_deriv( self, x ) result( y )

    class(sll_t_spline_1d_non_uniform), intent(in) :: self
    sll_real64            , intent(in) :: x
    sll_real64 :: y

    sll_int32  :: icell
    sll_int32  :: ib
    sll_int32  :: j
    sll_real64 :: values(self%deg+1)

    icell = sll_f_find_cell( self%bsp, x )
    call sll_s_bsplines_eval_deriv( self%bsp, icell, x, values )
    y = 0.0_f64
    do j = 1, self%deg+1
      ib = mod( icell+j-2-self%offset+self%n, self%n ) + 1
      y = y + values(j)*self%bcoef(ib)
    end do

  end function f_spline_1d_non_uniform__eval_deriv

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate spline S at given 1D array of points x(:)
  !> @param[in]  self spline object
  !> @param[in]  x    1D array of evaluation points
  !> @param[out] y    1D array of spline values y(:)=S(x(:))
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_1d_non_uniform__eval_array( self, x, y )

    class(sll_t_spline_1d_non_uniform), intent(in   ) :: self
    sll_real64            , intent(in   ) :: x(:)
    sll_real64            , intent(  out) :: y(size(x))

    sll_int32 :: i

    do i = 1, size(x)
      y(i) = f_spline_1d_non_uniform__eval( self, x(i) )
    end do

  end subroutine s_spline_1d_non_uniform__eval_array

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate spline derivative S' at given 1D array of points x(:)
  !> @param[in]  self spline object
  !> @param[in]  x    1D array of evaluation points
  !> @param[out] y    1D array of spline derivatives y(:)=S'(x(:))
  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_1d_non_uniform__eval_array_deriv( self, x, y )

    class(sll_t_spline_1d_non_uniform), intent(in   ) :: self
    sll_real64            , intent(in   ) :: x(:)
    sll_real64            , intent(  out) :: y(size(x))

    sll_int32 :: i

    do i = 1, size(x)
      y(i) = f_spline_1d_non_uniform__eval_deriv( self, x(i) )
    end do

  end subroutine s_spline_1d_non_uniform__eval_array_deriv

  !-----------------------------------------------------------------------------
  !> @brief Destructor. Frees all memory in object (pointers, allocatables)
  !> @param[inout] self object to be freed
  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_non_uniform__free( self )

    class(sll_t_spline_1d_non_uniform), intent(inout) :: self

    ! deallocate arrays
    deallocate( self%bcoef  )
    deallocate( self%tau    )
    if (allocated( self%q )) deallocate( self%q      )
    deallocate( self%bsdx   )

    ! free attribute objects
    call sll_s_bsplines_free( self%bsp )
    call schur_complement_free( self%schur )

  end subroutine s_spline_1d_non_uniform__free

end module sll_m_spline_1d_non_uniform
