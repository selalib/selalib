!> @ingroup splines
!> Implements arbitrary degree bspline interpolation on a uniform grid
!> given a B-Spline object from sll_m_arbitrary_degree_splines
module sll_m_bspline_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
     sll_p_periodic, &
     sll_p_hermite, &
     sll_p_greville, &
     sll_p_mirror

use sll_m_arbitrary_degree_splines, only: &
     sll_t_arbitrary_degree_spline_1d, &
     sll_s_arbitrary_degree_spline_1d_init, &
     sll_s_arbitrary_degree_spline_1d_free, &
     sll_f_find_cell, &
     sll_s_splines_at_x, &
     sll_s_spline_derivatives_at_x, &
     sll_s_splines_and_n_derivs_at_x, &
     sll_s_uniform_b_splines_at_x

use schur_complement, only: &
  schur_complement_solver, &
  schur_complement_fac   , &
  schur_complement_slv   , &
  schur_complement_free

implicit none

public :: &
  sll_t_bspline_1d,                     &
  sll_s_bspline_1d_init,                &
  sll_s_bspline_1d_free,                &
  sll_s_bspline_1d_compute_interpolant, &
  sll_f_bspline_1d_eval,                & ! scalar functions for evaluation
  sll_f_bspline_1d_eval_deriv,          &
  sll_s_bspline_1d_eval_array,          & ! vector subroutines for evaluation
  sll_s_bspline_1d_eval_array_deriv,    &
  sll_f_bspline_1d_get_coeff

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief
!> basic type for one-dimensional B-spline interpolation.
type :: sll_t_bspline_1d

  type (sll_t_arbitrary_degree_spline_1d) :: bsp
  sll_int32                               :: n        ! dimension of spline space
  sll_int32                               :: deg      ! degree of spline
  sll_real64, allocatable                 :: tau(:)   ! interpolation points
  sll_real64, allocatable                 :: bcoef(:) ! bspline coefficients
  sll_int32                               :: bc_type  ! boundary condition
  sll_real64, allocatable                 :: q(:,:)   ! triangular factorization
  sll_real64, allocatable                 :: values(:)
  sll_real64, allocatable                 :: bsdx(:,:)
  sll_int32,  allocatable                 :: ipiv(:)
  sll_real64                              :: length
  sll_int32                               :: offset
  type(schur_complement_solver)           :: schur
  sll_real64, allocatable                 :: bc_left (:)
  sll_real64, allocatable                 :: bc_right(:)

end type sll_t_bspline_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function sll_f_bspline_1d_get_coeff( self ) result( ptr )

    type(sll_t_bspline_1d), target, intent(in) :: self
    sll_real64, pointer :: ptr(:)

    ptr => self%bcoef

  end function sll_f_bspline_1d_get_coeff

  !> @brief Constructor for sll_t_bspline_1d object
  !> @param[in] num_pts  Number of points where the data to be
  !>                     interpolated are represented.
  !> @param[in] degree   Spline degree
  !> @param[in] xmin     Minimum value of the abscissae where data are meant
  !>                     to be interpolated.
  !> @param[in] xmax     Maximum value of the abscissae where data are meant
  !>                     to be interpolated.
  !> @param[in] bc_xmin  Boundary condition at x=xmin
  !> @param[in] bc_xmax  Boundary condition at x=xmax
  subroutine sll_s_bspline_1d_init( &
    self,           &
    num_pts,        &
    degree,         &
    xmin,           &
    xmax,           &
    bc_xmin,        &
    bc_xmax,        &
    ! optional arguments (TODO: verify if they should be removed)
    spline_bc_type, &
    bc_left,        &
    bc_right )

    type(sll_t_bspline_1d), intent(  out) :: self
    sll_int32             , intent(in   ) :: num_pts
    sll_int32             , intent(in   ) :: degree
    sll_real64            , intent(in   ) :: xmin
    sll_real64            , intent(in   ) :: xmax
    sll_int32             , intent(in   ) :: bc_xmin
    sll_int32             , intent(in   ) :: bc_xmax
    sll_int32 , optional  , intent(in   ) :: spline_bc_type
    sll_real64, optional  , intent(in   ) :: bc_left (:)
    sll_real64, optional  , intent(in   ) :: bc_right(:)

    sll_int32, parameter :: &
      bc_types_allowed(*) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

    sll_real64 :: grid(num_pts)
    sll_int32  :: i, j
    sll_real64 :: delta
    sll_int32  :: bc_type

    ! Sanity checks
    SLL_ASSERT( num_pts >= 2 )
    SLL_ASSERT( degree  >= 1 )
    SLL_ASSERT( xmin < xmax )
    SLL_ASSERT( any( bc_xmin == bc_types_allowed ) )
    SLL_ASSERT( any( bc_xmax == bc_types_allowed ) )

    ! NOTE: in the future different boundary conditions at xmin and xmax
    !       should be considered. For now we only check that bc_xmin==bc_xmax
    SLL_ASSERT( bc_xmin == bc_xmax )
    bc_type = bc_xmin

    ! knot point mirroring works only for Hermite boundary conditions. Check this
    if (present(spline_bc_type).and.(spline_bc_type == sll_p_mirror)) then
       if (.not.(bc_type == sll_p_hermite)) then
          print*, 'Mirror knot sequence at boundary only works with Hermite BC'
          stop
       end if
    end if

    ! set first attributes
    self%bc_type = bc_type
    if (bc_type == sll_p_periodic) then
       self%n = num_pts - 1 ! dimension of periodic spline space
       self%offset = degree/2  ! offset needed for periodic spline evaluation
    else
       self%n = num_pts + degree - 1 ! dimension of non periodic spline space
       self%offset = 0
    end if
    self%deg     = degree
    self%length  = xmax - xmin
    ! define grid for knots
    delta = self%length / real(num_pts-1,f64)
    do i=1,num_pts-1
       grid(i) = xmin + real(i-1,f64) * delta
    end do
    grid(num_pts) = xmax
    ! construct a sll_t_arbitrary_degree_spline_1d object
    if (present(spline_bc_type)) then
       call sll_s_arbitrary_degree_spline_1d_init(self%bsp, degree, grid,  &
            num_pts, bc_type, spline_bc_type )
    else
       call sll_s_arbitrary_degree_spline_1d_init(self%bsp, degree, grid,  &
            num_pts, bc_type)
    end if

    allocate( self%bcoef (self%n) )
    allocate( self%values(self%deg+1) )
    allocate( self%bsdx  (self%deg/2+1,self%deg+1) )
    allocate( self%ipiv  (self%n) )

    select case (bc_type)
    case (sll_p_periodic)
      ! Allocate tau which contains interpolation points
      allocate( self%tau(self%n) )
      if (modulo(degree,2) == 0) then
        ! for even degree interpolation points are cell midpoints
        do i = 1, self%n
          self%tau(i) = xmin + (real(i,f64)-0.5_f64) * delta
        end do
      else
        ! for odd degree interpolation points are grid points excluding last point
        do i = 1, self%n
          self%tau(i) = xmin + real(i-1,f64) * delta
        end do
      end if
    case (sll_p_hermite)
      ! In this case only interior points are saved in tau
       if (modulo(degree,2) == 0) then
         ! Allocate tau which contains interior interpolation points
         allocate( self%tau(num_pts-1) )
         ! for even degree interpolation points are cell midpoints
         do i = 1, num_pts - 1
           self%tau(i) = xmin + (real(i,f64)-0.5_f64) * delta
         end do
       else
         ! Allocate tau which contains interior interpolation points
         allocate( self%tau(num_pts) )
          ! for odd degree interpolation points are grid points including last point
          do i = 1, num_pts
             self%tau(i) = xmin + real(i-1,f64) * delta
          end do
       end if
       ! deal with possible boundary conditions
       if (present(bc_left).and. present(bc_right)) then
          SLL_ASSERT((size(bc_left )==self%deg/2))
          SLL_ASSERT((size(bc_right)==self%deg/2))
          allocate( self%bc_left (self%deg/2) )
          allocate( self%bc_right(self%deg/2) )
          self%bc_left  = bc_left
          self%bc_right = bc_right
       end if
    case (sll_p_greville)
      ! Allocate tau which contains interpolation points
      allocate( self%tau(self%n) )
      do i = 1, self%n
        ! Set interpolation points to be Greville points
        self%tau(i) = 0.0_f64
        do j=1-degree,0
          self%tau(i) = self%tau(i) + self%bsp%knots(i+j)
        end do
        self%tau(i) = self%tau(i) / real(degree,f64)
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

  end subroutine sll_s_bspline_1d_init

  !> @brief private subroutine for assembling and factorizing
  !> linear system needed for periodic spline interpolation
  !> @param[inout] self bspline interpolation object
  subroutine build_system_periodic( self )
    type(sll_t_bspline_1d), intent(inout) :: self

    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: sizeq

    ! evaluate bsplines at interpolation points
    if (modulo(self%deg,2) == 0) then
       ! for even degree interpolation points are cell midpoints
       call sll_s_uniform_b_splines_at_x(self%deg, 0.5_f64, self%values)
    else
       ! for odd degree interpolation points are grid points
       call sll_s_uniform_b_splines_at_x(self%deg, 0.0_f64, self%values)
    end if
    sizeq = self%deg/2
    ! assemble matrix q for linear system
    allocate( self%q(2*sizeq+1,self%n) )

    do j=1,self%n
       do i=1, 2*sizeq+1
          self%q(i,j) = self%values(i)
       end do
    end do

    ! Perform LU decomposition of matrix q
    call schur_complement_fac( self%schur, self%n, sizeq, self%q )

  end subroutine build_system_periodic

  !> @brief private subroutine for assembling and factorizing
  !> linear system needed for spline interpolation at Greville points
  !> @param[inout] self bspline interpolation object
  subroutine build_system_greville( self )
    type(sll_t_bspline_1d), intent(inout) :: self

    sll_int32    :: iflag
    sll_int32    :: j
    sll_int32    :: k
    sll_int32    :: ii
    sll_int32    :: jj
    sll_int32    :: icell
    sll_real64   :: x

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
    call sll_s_splines_at_x(self%bsp, icell, x, self%values)
    ! iterate only to k+1 as last bspline is 0
    do j=1,k+1
       jj = icell+j-1
       self%q(ii-jj+k+1,jj) = self%values(j)
    end do

    do ii=2,self%n - 1
       x = self%tau(ii)
       icell = sll_f_find_cell(self%bsp, x )
       call sll_s_splines_at_x(self%bsp, icell, x, self%values)
       do j=1,k+2
          jj = icell+j-1
          self%q(ii-jj+k+1,jj) = self%values(j)
       end do
    end do
    ! Treat i=self%n separately
    x =  self%bsp%xmax
    icell = self%n - self%deg
    ii = self%n  ! last Greville point
    call sll_s_splines_at_x(self%bsp, icell, x, self%values)
    ! iterate only to k as first bspline is 0
    do j=2,k+2
       jj = icell+j-1
       self%q(ii-jj+k+1,jj) = self%values(j)
    end do
    ! Perform LU decomposition of matrix q
    call banfac ( self%q, 2*k+1, self%n, k, k, iflag )

  end subroutine build_system_greville

  !> @brief private subroutine for assembling and factorizing
  !> linear system needed for spline interpolation with Hermite boundary conditions
  !> @param[inout] self bspline interpolation object
  !> @param[in] nder number of derivatives to be computed
  subroutine build_system_with_derivative( self )

    type(sll_t_bspline_1d), intent(inout) :: self

    sll_int32   :: nbc
    sll_int32   :: iflag
    sll_int32   :: i
    sll_int32   :: j
    sll_int32   :: ii
    sll_int32   :: jj
    sll_int32   :: offset
    sll_int32   :: k
    sll_int32   :: icell
    sll_real64  :: x

    ! number of boundary conditions needed depending on spline degree
    k = self%deg-1
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
    call sll_s_splines_and_n_derivs_at_x( self%bsp, icell, x , nbc , self%bsdx)
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
       icell = sll_f_find_cell(self%bsp, x )
       call sll_s_splines_at_x(self%bsp, icell, x, self%values)
       do j=1,k+2
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = self%values(j)
       end do
    end do
    ! boundary conditions at xmax
    x = self%bsp%xmax
    icell = self%n - self%deg
    call sll_s_splines_and_n_derivs_at_x( self%bsp, icell, x , nbc , self%bsdx)
    do i= 1, nbc
       ii = ii + 1
       do j=2,k+2
          jj = icell+j-1
          self%q( ii-jj+2*k+1,jj) = self%bsdx(i+offset,j)
       end do
    end do

    ! Perform LU decomposition of matrix q with Lapack
    call dgbtrf(self%n,self%n,k,k,self%q,3*k+1,self%ipiv,iflag)

  end subroutine build_system_with_derivative

  !> @brief subroutine for computing interpolating spline
  !> Computes 1D spline interpolating function gtau at interpolation points
  !> @param[inout] self bspline interpolation object
  !> @param[in] val_min (optional) array containing boundary conditions at x1_min
  !> @param[in] val_max (optional) array containing boundary conditions at x1_max
  !>
  subroutine sll_s_bspline_1d_compute_interpolant( self, gtau, val_min, val_max )

    type(sll_t_bspline_1d), intent(inout) :: self
    sll_real64            , intent(in   ) :: gtau(:)
    sll_real64, optional  , intent(in   ) :: val_min(:)
    sll_real64, optional  , intent(in   ) :: val_max(:)

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
       call schur_complement_slv ( self%schur, self%n, self%deg/2, self%q,  self%bcoef )
    case (sll_p_greville)
       self%bcoef = gtau
       call banslv ( self%q, 2*self%deg-1, self%n, self%deg-1, self%deg-1, self%bcoef )
    case (sll_p_hermite)
       ! number of needed conditions at boundary
       ncond = self%deg/2
       if (present(val_min)) then
          self%bcoef(1:ncond) = val_min(1:ncond)
       else  ! set needed boundary values to 0
          self%bcoef(1:ncond) = 0.0_f64
       end if
       self%bcoef(ncond+1:self%n-ncond) = gtau(1:self%n-2*ncond)
       if (present(val_max)) then
          self%bcoef(self%n-ncond+1:self%n) = val_max(1:ncond)
       else ! set needed boundary values to 0
          self%bcoef(self%n-ncond+1:self%n) = 0.0_f64
       end if
       ! Use Lapack to solve banded system
       k = self%deg-1
       call dgbtrs('N',self%n,k,k,1,self%q,3*k+1,self%ipiv,self%bcoef, self%n, iflag)
    end select

  end subroutine sll_s_bspline_1d_compute_interpolant

!> @brief returns the value of the image of an abscissae,
!> The spline coefficients
!> used are stored in the spline object pointer.
!> @param[in] x input double-precison element containing the
!> abscissae to be interpolated.
!> @param[out] y output double-precision element containing the
!> results of the interpolation.
!> @param[inout] spline the spline object pointer, duly initialized and
!> already operated on by the sll_s_bspline_1d_compute_interpolant() subroutine.
  function sll_f_bspline_1d_eval( self, x ) result( y )

    type(sll_t_bspline_1d), intent(in) :: self
    sll_real64            , intent(in) :: x
    sll_real64 :: y

    sll_int32 :: icell
    sll_int32 :: ib
    sll_int32 :: j

    ! get bspline values at x
    icell =  sll_f_find_cell( self%bsp, x )
    call sll_s_splines_at_x(self%bsp, icell, x, self%values)
    y = 0.0_f64
    do j=1,self%deg+1
       ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
       y = y + self%values(j)*self%bcoef(ib)
    enddo

  end function sll_f_bspline_1d_eval

!> @brief returns the values of the images of a collection of
!> abscissae, represented by a 1D array in another output array.
!> The spline coefficients used are stored in the spline object pointer.
!> @param[in] x input double-precison element array containing the
!> abscissae to be interpolated.
!> @param[out] y output double-precision element array containing the
!> results of the interpolation.
!> @param[in] n the number of elements of the input array which are to be
!> interpolated.
!> @param[inout] spline the spline object pointer, duly initialized and
!> already operated on by the sll_s_bspline_1d_compute_interpolant() subroutine.
subroutine sll_s_bspline_1d_eval_array( self, n, x, y )

  type(sll_t_bspline_1d), intent(in   ) :: self
  sll_int32             , intent(in   ) :: n
  sll_real64            , intent(in   ) :: x(n)
  sll_real64            , intent(  out) :: y(n)

  sll_real64 :: val
  sll_real64 :: xx
  sll_int32  :: i
  sll_int32  :: j
  sll_int32  :: ib
  sll_int32  :: icell

  do i= 1, n
     xx = x(i)
     val = 0.0_f64
     ! get bspline values at x
     icell =  sll_f_find_cell( self%bsp, xx )
     call sll_s_splines_at_x(self%bsp, icell, xx, self%values)
     do j=1, self%deg+1
        ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
        val = val + self%values(j)*self%bcoef(ib)
     enddo
     y(i) = val
  end do

end subroutine sll_s_bspline_1d_eval_array

!> @brief returns the values of the derivatives evaluated at a point.
!> The spline coefficients used are stored in the spline object.
!> @param[inout] self the spline object pointer, duly initialized and
!> already operated on by the sll_s_bspline_1d_compute_interpolant() subroutine.
!> @param[in] x  abscissa to be interpolated.
!> @param[out] y  result of the interpolation.
function sll_f_bspline_1d_eval_deriv( self, x ) result( y )

  type(sll_t_bspline_1d), intent(in) :: self
  sll_real64            , intent(in) :: x
  sll_real64 :: y

  sll_int32 :: icell
  sll_int32 :: ib
  sll_int32 :: j

  ! get bspline derivatives at x
  icell =  sll_f_find_cell( self%bsp, x )
  call sll_s_spline_derivatives_at_x(self%bsp, icell, x, self%values)
  y = 0.0_f64
  do j=1,self%deg+1
     ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
     y = y + self%values(j)*self%bcoef(ib)
  enddo

end function sll_f_bspline_1d_eval_deriv

!> @brief returns the values of the derivatives evaluated at a
!> collection of abscissae stored by a 1D array in another output
!> array. The spline coefficients used are stored in the spline
!> object pointer.
!> @param[in] x input double-precison element array containing the
!> abscissae to be interpolated.
!> @param[out] y output double-precision element array containing the
!> results of the interpolation.
!> @param[in] n the number of elements of the input array which are to be
!> interpolated.
!> @param[inout] spline the spline object pointer, duly initialized and
!> already operated on by the sll_s_bspline_1d_compute_interpolant() subroutine.
subroutine sll_s_bspline_1d_eval_array_deriv( self, n, x, y )

  type(sll_t_bspline_1d), intent(in   ) :: self
  sll_int32             , intent(in   ) :: n
  sll_real64            , intent(in   ) :: x(n)
  sll_real64            , intent(  out) :: y(n)

  sll_real64 :: val
  sll_real64 :: xx
  sll_int32  :: i
  sll_int32  :: j
  sll_int32  :: ib
  sll_int32  :: icell

  do i= 1, n
     xx = x(i)
     val = 0.0_f64
     ! get bspline derivatives at xx
     icell =  sll_f_find_cell( self%bsp, xx )
     call sll_s_spline_derivatives_at_x(self%bsp, icell, xx, self%values)
     do j=1, self%deg+1
        ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
        val = val + self%values(j)*self%bcoef(ib)
     enddo
     y(i) = val
  end do

end subroutine sll_s_bspline_1d_eval_array_deriv

!> @brief Destructor. Frees all pointers of object
!> @param[inout] self object to be freed
subroutine sll_s_bspline_1d_free( self )

  type(sll_t_bspline_1d), intent(inout) :: self

  ! deallocate arrays
  deallocate( self%bcoef  )
  deallocate( self%tau    )
  deallocate( self%q      )
  deallocate( self%values )
  deallocate( self%bsdx   )

  ! free attribute objects
  call sll_s_arbitrary_degree_spline_1d_free( self%bsp )
  call schur_complement_free( self%schur )

end subroutine sll_s_bspline_1d_free

end module sll_m_bspline_1d
