!> @ingroup splines
!> @brief   Arbitrary degree spline interpolation on uniform grid
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_1d_uniform
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite , &
    sll_p_greville

  use sll_m_spline_1d_base, only: &
    sll_c_spline_1d

  use sll_m_bsplines, only: &
    sll_s_uniform_bsplines_eval_basis, &
    sll_s_uniform_bsplines_eval_deriv, &
    sll_s_uniform_bsplines_eval_basis_and_n_derivs

  use schur_complement, only: &
    schur_complement_solver, &
    schur_complement_fac   , &
    schur_complement_slv   , &
    schur_complement_free

  implicit none

  public :: &
    sll_t_spline_1d_uniform

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Allowed boundary conditions
  integer, parameter :: &
    allowed_bcs(*) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  !> 1D spline on uniform grid
  type, extends(sll_c_spline_1d) :: sll_t_spline_1d_uniform

    integer , private :: deg     ! spline degree (= max order of piecewise polynomial)
    integer , private :: n       ! dimension of spline space

    real(wp), private :: xmin
    real(wp), private :: xmax
    integer , private :: ncells
    real(wp), private :: dx

    integer , private :: bc_xmin ! boundary condition type at x=xmin
    integer , private :: bc_xmax ! boundary condition type at x=xmax

    integer , private :: offset  ! shift in bsplines indexing (periodic only)

    real(wp), private, allocatable :: bcoef(:) ! spline coefficients (w.r.t. basis)
    real(wp), private, allocatable :: tau(:)   ! interpolation points

    ! Local storage for linear system in interpolation problem
    real(wp),         allocatable, private :: q(:,:)
    integer ,         allocatable, private :: ipiv(:)
    type(schur_complement_solver), private :: schur

  contains

    procedure :: init                => s_spline_1d_uniform__init
    procedure :: free                => s_spline_1d_uniform__free
    procedure :: compute_interpolant => s_spline_1d_uniform__compute_interpolant
    procedure :: eval                => f_spline_1d_uniform__eval
    procedure :: eval_deriv          => f_spline_1d_uniform__eval_deriv
    procedure :: eval_array          => s_spline_1d_uniform__eval_array
    procedure :: eval_array_deriv    => s_spline_1d_uniform__eval_array_deriv
    procedure :: get_coeff           => f_spline_1d_uniform__get_coeff
    procedure :: get_interp_points   => s_spline_1d_uniform__get_interp_points

  end type sll_t_spline_1d_uniform

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_get_icell_and_offset( xmin, xmax, ncells, x, icell, offset )

    real(wp), intent(in   ) :: xmin
    real(wp), intent(in   ) :: xmax
    integer , intent(in   ) :: ncells
    real(wp), intent(in   ) :: x
    integer , intent(  out) :: icell
    real(wp), intent(  out) :: offset

    real(wp) :: x_normalized  ! 0 <= x_normalized <= num_cells

    SLL_ASSERT( x >= xmin )
    SLL_ASSERT( x <= xmax )

    if (x == xmin) then
      icell  = 1
      offset = 0.0_wp
    else if (x == xmax) then
      icell  = ncells
      offset = 1.0_wp
    else
      x_normalized = (x-xmin) / (xmax-xmin) * real(ncells,wp)
      icell        = int( x_normalized )
      offset       = x_normalized - real( icell, wp )
      icell        = icell + 1
    end if

  end subroutine s_get_icell_and_offset

  !-----------------------------------------------------------------------------
  function f_spline_1d_uniform__get_coeff( self ) result( ptr )

    class(sll_t_spline_1d_uniform), target, intent(in) :: self
    real(wp), pointer :: ptr(:)

    ptr => self%bcoef

  end function f_spline_1d_uniform__get_coeff

  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_uniform__get_interp_points( self, tau )

    class(sll_t_spline_1d_uniform), intent(in   ) :: self
    real(wp),          allocatable, intent(  out) :: tau(:)

    SLL_ASSERT( allocated( self%tau ) )
    allocate( tau(size(self%tau)), source=self%tau )

  end subroutine s_spline_1d_uniform__get_interp_points

  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_uniform__init( &
      self   , &
      degree , &
      ncells , &
      xmin   , &
      xmax   , &
      bc_xmin, &
      bc_xmax )

    class(sll_t_spline_1d_uniform), intent(  out) :: self
    integer                       , intent(in   ) :: degree
    integer                       , intent(in   ) :: ncells
    real(wp)                      , intent(in   ) :: xmin
    real(wp)                      , intent(in   ) :: xmax
    integer                       , intent(in   ) :: bc_xmin
    integer                       , intent(in   ) :: bc_xmax

    character(len=*), parameter :: this_sub_name = "sll_t_spline_1d_uniform % init"
    character(len=256) :: err_msg

    real(wp) :: shift
    integer  :: i
    integer  :: size_tau
    real(wp), allocatable :: knots(:)

    ! Sanity checks
    SLL_ASSERT( degree >= 1 )
    SLL_ASSERT( ncells >= 1 )
    SLL_ASSERT( xmin < xmax )
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )

    if (any( [bc_xmin,bc_xmax]== sll_p_periodic ) .and. bc_xmin /= bc_xmax) then
      err_msg = "Periodic boundary conditions mismatch: bc_xmin /= bc_xmax."
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if

    ! TODO: implement different BCs at xmin/xmax
    if (bc_xmin /= bc_xmax) then
      err_msg = "Different BCs at xmin/xmax not yet implemented."
      SLL_ERROR( this_sub_name, trim( err_msg ) )
    end if

    ! Compute cell size in uniform grid
    self%dx = (xmax-xmin) / real(ncells,wp)

    ! Determine number of degrees of freedom
    ! Determine offset (non-zero for periodic spline evaluation)
    if (bc_xmin == sll_p_periodic) then
      self%n      = ncells
      self%offset = degree/2
    else
      self%n      = ncells+degree
      self%offset = 0
    end if

    ! Store other info
    self%deg     = degree
    self%ncells  = ncells
    self%xmin    = xmin
    self%xmax    = xmax
    self%bc_xmin = bc_xmin
    self%bc_xmax = bc_xmax

    ! Allocate array of spline coefficients
    allocate( self%bcoef(self%n) )

    ! Determine interpolation points
    ! TODO: handle different BCs at xmin/xmax
    select case (bc_xmin)
 
    ! 1) Periodic case: 
    !    - for even degree use cell midpoints
    !    - for odd  degree use grid points excluding last point
    case (sll_p_periodic)

      if (modulo(degree,2) == 0) then
        shift    = 0.5_wp
        size_tau = self%n
      else
        shift    = 0.0_wp
        size_tau = self%n
      end if
      allocate( self%tau (1:size_tau) )
      self%tau = [(xmin + (real(i-1,wp)+shift)*self%dx, i=1,size_tau)]

    ! 2) Hermite-Hermite case:
    !    - for even degree use cell midpoints
    !    - for odd  degree use all grid points (including last point)
    case (sll_p_hermite)

       if (modulo(degree,2) == 0) then
         shift    = 0.5_wp
         size_tau = ncells
       else
         shift    = 0.0_wp
         size_tau = ncells+1
       end if
       allocate( self%tau (1:size_tau) )
       self%tau = [(xmin + (real(i-1,wp)+shift)*self%dx, i=1,size_tau)]

    ! 3) Greville-Greville case:
    !    - for even degree use cell midpoints
    !    - for odd  degree use all grid points (including last point)
    !    - calculate additional interpolation points in first and last cells
    case (sll_p_greville)

      allocate( self%tau(self%n) )

      associate( num_pts => ncells+1 )
        allocate( knots (1-degree:num_pts+degree) )
        knots(1-degree:0) = xmin
        knots(1:num_pts ) = [(xmin+real(i-1,wp)*self%dx, i=1,num_pts)]
        knots(num_pts+1:num_pts+degree) = xmax
        do i = 1, self%n
          self%tau(i) = sum( knots(i+1-degree:i) ) / real(degree,wp)
        end do
        deallocate( knots )
      end associate

    end select

    ! Special case: linear spline
    ! No need for matrix assembly
    if (self%deg == 1) return

    ! Assemble matrix $M_{ij} = B_j(\tau_i)$ for spline interpolation
    ! TODO: handle mixed Hermite/Greville BCs
    select case (bc_xmin)
    case (sll_p_periodic); call s_build_system_periodic( self )
    case (sll_p_hermite ); call s_build_system_hermite ( self )
    case (sll_p_greville); call s_build_system_greville( self )
    end select

  end subroutine s_spline_1d_uniform__init

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system needed for periodic spline interpolation
  !> @param[inout] self spline 1D object
  !-----------------------------------------------------------------------------
  subroutine s_build_system_periodic( self )

    class(sll_t_spline_1d_uniform), intent(inout) :: self

    integer  :: j
    integer  :: k
    real(wp) :: offset
    real(wp) :: values(self%deg+1)

    ! Evaluate B-splines at interpolation points, noting that:
    !   a) points are regularly spaced
    !   b) basis is translationally invariant
    ! For even degree interpolation points are cell midpoints
    ! For odd  degree interpolation points are grid points
    if (modulo(self%deg,2) == 0) then
      offset = 0.5_wp
    else
      offset = 0.0_wp
    end if
    call sll_s_uniform_bsplines_eval_basis( self%deg, offset, values )

    ! Allocate matrix Q for linear system:
    !   - M is circulant banded matrix (sparse)
    !   - Q is compressed form of A (circulant diagonals of M are rows of Q)
    k = self%deg/2
    allocate( self%q(2*k+1,self%n) )

    ! Assemble matrix Q:
    !  - M has constant elements on each diagonal
    !  - Q has constant elements on each row
    do j = 1, self%n
      self%q(:,j) = values(1:2*k+1)
    end do

    ! Prepare solution of linear system by performing:
    ! a) Block decomposition of M = [[A,B],[C,D]] with
    !     - A (n-k)x(n-k) banded Toeplitz matrix (square, large)
    !     - B (n-k)x(  k) low-rank Toeplitz (rectangular, empty in the center)
    !     - C (  k)x(n-k) low-rank Toeplitz (rectangular, empty in the center)
    !     - D (  k)x(  k) fully dense Toeplitz matrix (square, small)
    ! b) LU factorization of matrix A
    ! c) Calculation of Schur complement of A: H = D - C A^{-1} B
    ! d) LU factorization of Schur complement H
    call schur_complement_fac( self%schur, self%n, k, self%q )

  end subroutine s_build_system_periodic

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for spline interpolation with Hermite BCs
  !> @param[inout] self spline 1D object
  !-----------------------------------------------------------------------------
  subroutine s_build_system_hermite( self )

    class(sll_t_spline_1d_uniform), intent(inout) :: self

    integer  :: nbc
    integer  :: iflag
    integer  :: i
    integer  :: j
    integer  :: ii
    integer  :: jj
    integer  :: offset
    integer  :: k
    integer  :: icell
    real(wp) :: x
    real(wp) :: x_offset

    real(wp) :: values(self%deg+1)
    real(wp) :: derivs(0:self%deg/2,1:self%deg+1)


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
    self%q = 0.0_wp

    ! For even degree splines interpolation points are at cell midpoints
    ! value of function on the boundary needed as additional boundary conditions
    ! for odd degree splines  baundary is in the interpolation points
    ! only derivative values are needed as boundary conditions
    if (modulo( self%deg, 2 ) == 0) then ! spline degree even
       offset = 1
    else
       offset = 0
    end if

    ! boundary conditions at xmin
    ii       = 0
    x        = self%xmin
    icell    = 1
    x_offset = 0.0_wp
    call sll_s_uniform_bsplines_eval_basis_and_n_derivs( self%deg, x_offset, nbc, derivs )
    ! When using uniform B-splines, the cell size is normalized between 0 and 1,
    ! hence the i-th derivative must be divided by dx^i to scale back the result
    do i = 1, nbc-offset
      derivs(i,:) = derivs(i,:) / self%dx**i
    end do
    do i = 1, nbc
       ! iterate only to k+1 as last bspline is 0
       ii = ii + 1
       do j = 1, k+1
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = derivs(i-offset,j)
       end do
    end do

    ! interpolation points
    do i = 1, self%n-2*nbc
       ii = ii + 1
       x = self%tau(i)
       call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, x_offset )
       call sll_s_uniform_bsplines_eval_basis( self%deg, x_offset, values )
       do j=1,k+2
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = values(j)
       end do
    end do

    ! boundary conditions at xmax
    x        = self%xmax
    icell    = self%n-self%deg
    x_offset = 1.0_wp
    call sll_s_uniform_bsplines_eval_basis_and_n_derivs( self%deg, x_offset, nbc, derivs )
    ! When using uniform B-splines, the cell size is normalized between 0 and 1,
    ! hence the i-th derivative must be divided by dx^i to scale back the result
    do i = 1, nbc-offset
      derivs(i,:) = derivs(i,:) / self%dx**i
    end do
    do i = 1, nbc
       ii = ii + 1
       do j = 2, k+2
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = derivs(i-offset,j)
       end do
    end do

    ! Perform LU decomposition of matrix q with Lapack
    allocate( self%ipiv (self%n) )
    call dgbtrf( self%n, self%n, k, k, self%q, 3*k+1, self%ipiv, iflag )

  end subroutine s_build_system_hermite

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for spline interpolation at Greville points
  !> @param[inout] self bspline object
  !-----------------------------------------------------------------------------
  subroutine s_build_system_greville( self )

    class(sll_t_spline_1d_uniform), intent(inout) :: self

    integer  :: iflag
    integer  :: j
    integer  :: k
    integer  :: ii
    integer  :: jj
    integer  :: icell
    real(wp) :: x
    real(wp) :: offset
    real(wp) :: values(self%deg+1)

    ! The matrix q is a banded matrix using the storage required by banfac (De Boor)
    ! It has k bands above diagonal, k bands below and the diagonal itself
    ! The term A(ii,jj) of the full matrix is stored in q(ii-jj+k+1,jj)
    ! The Bspline interpolation matrix at Greville points x_ii is
    ! A(ii,jj) = B_jj(x_ii)
    k = self%deg - 1
    allocate( self%q(2*k+1,self%n) )
    self%q = 0.0_f64

    ! Treat i=1 separately
    x     =  self%xmin
    icell = 1
    ii    = 1   ! first Greville point
    call sll_s_uniform_bsplines_eval_basis( self%deg, 0.0_wp, values )
    ! iterate only to k+1 as last bspline is 0
    do j=1,k+1
       jj = icell+j-1
       self%q(ii-jj+k+1,jj) = values(j)
    end do

    do ii=2,self%n - 1
       x = self%tau(ii)
       call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, offset )
       call sll_s_uniform_bsplines_eval_basis( self%deg, offset, values )
       do j=1,k+2
          jj = icell+j-1
          self%q(ii-jj+k+1,jj) = values(j)
       end do
    end do

    ! Treat i=self%n separately
    x     = self%xmax
    icell = self%n - self%deg
    ii    = self%n  ! last Greville point
    call sll_s_uniform_bsplines_eval_basis( self%deg, 1.0_wp, values )
    ! iterate only to k as first bspline is 0
    do j=2,k+2
       jj = icell+j-1
       self%q(ii-jj+k+1,jj) = values(j)
    end do

    ! Perform LU decomposition of matrix q
    call banfac ( self%q, 2*k+1, self%n, k, k, iflag )

  end subroutine s_build_system_greville

  !-----------------------------------------------------------------------------
  !> @brief        Compute interpolating 1D spline
  !> Computes coefficients of 1D bspline that interpolates function on grid.
  !> If Hermite BCs are used, derivatives at boundary are also needed.
  !>
  !> @param[inout] self        1D bspline object
  !> @param[in]    gtau        function values of interpolation points
  !> @param[in]    derivs_xmin (optional) array with boundary conditions at xmin
  !> @param[in]    derivs_xmax (optional) array with boundary conditions at xmax
  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_uniform__compute_interpolant( self, gtau, derivs_xmin, derivs_xmax )

    class(sll_t_spline_1d_uniform), intent(inout) :: self
    real(wp)                      , intent(in   ) :: gtau(:)
    real(wp),             optional, intent(in   ) :: derivs_xmin(:)
    real(wp),             optional, intent(in   ) :: derivs_xmax(:)

    integer :: ncond, k, iflag

    ! Special case: linear spline
    if (self%deg == 1) then
      self%bcoef(:) = gtau(1:self%n)
      return
    end if

    ! Degree > 1: boundary conditions matter
    ! TODO: handle different BCs at xmin/xmax
    select case (self%bc_xmin)

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

  end subroutine s_spline_1d_uniform__compute_interpolant

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d_uniform__eval( self, x ) result( y )

    class(sll_t_spline_1d_uniform), intent(in) :: self
    real(wp)                      , intent(in) :: x
    real(wp) :: y

    integer  :: icell
    real(wp) :: offset
    real(wp) :: values(self%deg+1)

    ! TODO: save larger array of coefficients in case of periodic BCs
    !       to avoid loop and usage of "mod( , )" function!!
    integer  :: j, ib

    call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, offset )
    call sll_s_uniform_bsplines_eval_basis( self%deg, offset, values )

    ! FIXME
!    y = dot_product( self%bcoef(icell-self%deg:icell), values )
    y = 0.0_wp
    do j = 1, self%deg+1
       ib = mod( icell+j-2-self%offset+self%n, self%n ) + 1
       y = y + values(j)*self%bcoef(ib)
    end do

  end function f_spline_1d_uniform__eval

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d_uniform__eval_deriv( self, x ) result( y )

    class(sll_t_spline_1d_uniform), intent(in) :: self
    real(wp)                      , intent(in) :: x
    real(wp) :: y

    integer  :: icell
    real(wp) :: offset
    real(wp) :: values(self%deg+1)

    ! TODO: save larger array of coefficients in case of periodic BCs
    !       to avoid loop and usage of "mod( , )" function!!
    integer  :: j, ib

    call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, offset )
    call sll_s_uniform_bsplines_eval_deriv( self%deg, offset, values )

    ! FIXME
!    y = dot_product( self%bcoef(icell-self%deg:icell), values )
    y = 0.0_wp
    do j = 1, self%deg+1
       ib = mod( icell+j-2-self%offset+self%n, self%n ) + 1
       y = y + values(j)*self%bcoef(ib)
    end do

    ! When using uniform B-splines, the cell size is normalized between 0 and 1,
    ! hence the derivative must be divided by dx to scale back the result.
    y = y / self%dx

  end function f_spline_1d_uniform__eval_deriv

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_1d_uniform__eval_array( self, x, y )

    class(sll_t_spline_1d_uniform), intent(in   ) :: self
    real(wp)                      , intent(in   ) :: x(:)
    real(wp)                      , intent(  out) :: y(size(x))

    integer :: i

    do i = 1, size(x)
      y(i) = f_spline_1d_uniform__eval( self, x(i) )
    end do

  end subroutine s_spline_1d_uniform__eval_array

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_1d_uniform__eval_array_deriv( self, x, y )

    class(sll_t_spline_1d_uniform), intent(in   ) :: self
    real(wp)                      , intent(in   ) :: x(:)
    real(wp)                      , intent(  out) :: y(size(x))

    integer :: i

    do i = 1, size(x)
      y(i) = f_spline_1d_uniform__eval_deriv( self, x(i) )
    end do

  end subroutine s_spline_1d_uniform__eval_array_deriv

  !-----------------------------------------------------------------------------
  !> @brief Destructor. Frees all memory in object (pointers, allocatables)
  !> @param[inout] self object to be freed
  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_uniform__free( self )

    class(sll_t_spline_1d_uniform), intent(inout) :: self

    ! Deallocate arrays
    deallocate( self%bcoef  )
    deallocate( self%tau    )

    ! Arrays used by linear solvers
    if (allocated( self%q    )) deallocate( self%q    )
    if (allocated( self%ipiv )) deallocate( self%ipiv )

    ! Free attribute objects
    call schur_complement_free( self%schur )

  end subroutine s_spline_1d_uniform__free

end module sll_m_spline_1d_uniform
