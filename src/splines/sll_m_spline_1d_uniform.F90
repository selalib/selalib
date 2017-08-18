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

  use sll_m_spline_matrix, only: &
    sll_c_spline_matrix, &
    sll_s_spline_matrix_new

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

  !> 1D spline interpolation on uniform grid
  type, extends(sll_c_spline_1d) :: sll_t_spline_1d_uniform

    integer :: deg      ! spline degree (= order of piecewise polynomial)
    integer :: n        ! dimension of spline space
    integer :: bc_xmin  ! boundary condition type at x=xmin
    integer :: bc_xmax  ! boundary condition type at x=xmax
    integer :: offset   ! shift in B-spline indexing (periodic only)

    real(wp), private :: xmin     ! left  boundary coordinate
    real(wp), private :: xmax     ! right boundary coordinate
    integer , private :: nbc_xmin ! number of boundary conditions (derivatives) at x=xmin
    integer , private :: nbc_xmax ! number of boundary conditions (derivatives) at x=xmax
    integer , private :: ncells   ! number of cells
    logical , private :: even     ! true if deg even, false if deg odd
    integer , private :: mod      ! result of modulo(deg,2): 0 if deg even, 1 if deg odd

    real(wp), private :: dx       ! cell size
    real(wp), private :: inv_dx   ! inverse of cell size (=1/dx)

    real(wp), allocatable, private :: bcoef(:) ! B-splines' coefficients
    real(wp), allocatable, private :: tau(:)   ! interpolation points

    ! Polymorphic matrix object to store and solve linear system for interpolation
    class(sll_c_spline_matrix), allocatable, private :: matrix

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

    integer              :: i
    integer              :: ntau
    integer              :: isum
    integer, allocatable :: iknots(:)
    integer              :: kl, ku
    character(len=32)    :: matrix_type

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

    ! Compute cell size in uniform grid
    self%dx     = (xmax-xmin) / real(ncells,wp)
    self%inv_dx = 1.0_wp / self%dx

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
    self%deg      = degree
    self%ncells   = ncells
    self%xmin     = xmin
    self%xmax     = xmax
    self%bc_xmin  = bc_xmin
    self%bc_xmax  = bc_xmax
    self%nbc_xmin = merge( degree/2, 0, bc_xmin == sll_p_hermite )
    self%nbc_xmax = merge( degree/2, 0, bc_xmax == sll_p_hermite )
    self%mod      = modulo( degree, 2 )
    self%even     = (self%mod == 0)

    ! Allocate array of spline coefficients
    ! in case of periodic BCs, a larger array of coefficients is used in order
    ! to avoid a loop with calls to the "mod( , )" function at evaluation.
    if (self%bc_xmin == sll_p_periodic) then
      allocate( self%bcoef (1-degree:self%n+degree) )
    else
      allocate( self%bcoef (self%n) )
    end if

    !---------------------------------------------------------------------------
    ! Determine array tau of interpolation points
    !---------------------------------------------------------------------------

    ! Determine size of tau and allocate tau
    ntau = self%n - self%nbc_xmin - self%nbc_xmax
    allocate( self%tau(1:ntau) )

    if ( self%bc_xmin == sll_p_periodic ) then

      ! Periodic case:
      !   . for odd  degree, interpolation points are breakpoints (last excluded)
      !   . for even degree, interpolation points are cell centers
      associate( shift => merge( 0.5_wp, 0.0_wp, self%even ) )
        self%tau = [(xmin + (real(i-1,wp)+shift)*self%dx, i=1,ntau)]
      end associate

    else

      ! Non-periodic case: create array of temporary knots (integer shifts only)
      ! in order to compute interpolation points using Greville-style averaging:
      ! tau(i) = xmin + average(knots(i+1-degree:i)) * dx
      allocate( iknots (2-degree:ntau) )

      ! Additional knots near x=xmin
      associate( r => 2-degree, s => -self%nbc_xmin )
        select case (bc_xmin)
        case (sll_p_greville); iknots(r:s) = 0
        case (sll_p_hermite ); iknots(r:s) = [(i,i=r-s-1,-1)]
        end select
      end associate

      ! Knots inside the domain
      associate( r => -self%nbc_xmin+1, s => -self%nbc_xmin+1+self%ncells )
        iknots(r:s) = [(i,i=0,self%ncells)]
      end associate

      ! Additional knots near x=xmax
      associate( r => -self%nbc_xmin+1+self%ncells+1, s => ntau, n => self%ncells )
        select case (bc_xmax)
        case (sll_p_greville); iknots(r:s) = n
        case (sll_p_hermite ); iknots(r:s) = [(i,i=n+1,n+1+s-r)]
        end select
      end associate

      ! Compute interpolation points using Greville-style averaging
      associate( inv_deg => 1.0_wp / real( degree, wp ) )
        do i = 1, ntau
          isum = sum( iknots(i+1-degree:i) )
          if (modulo( isum, degree ) == 0) then
            self%tau(i) = xmin + real(isum/degree,wp) * self%dx
          else
            self%tau(i) = xmin + real(isum,wp) * inv_deg * self%dx
          end if
        end do
      end associate

      deallocate( iknots )

    end if

    !---------------------------------------------------------------------------
    ! Assemble dense matrix (B_j(tau(i))) for spline interpolation
    !---------------------------------------------------------------------------

    ! Special case: linear spline
    ! No need for matrix assembly
    if (self%deg == 1) return

    ! FIXME: In Hermite case ku and kl computed in general case when derivatives
    !        of B-splines do not vanish at boundary
    ! TODO: reduce kl and ku to take advantage of uniform grid
    select case( self%bc_xmin )
      case ( sll_p_periodic ); ku = (self%deg+1)/2
      case ( sll_p_hermite  ); ku = max( (self%deg+1)/2, self%deg-1 )
      case ( sll_p_greville ); ku = self%deg
    end select

    select case( self%bc_xmax )
      case ( sll_p_periodic ); kl = (self%deg+1)/2
      case ( sll_p_hermite  ); kl = max( (self%deg+1)/2, self%deg-1 )
      case ( sll_p_greville ); kl = self%deg
    end select

    if (self%bc_xmin == sll_p_periodic) then
      if (kl+1+ku >= self%n) then
        matrix_type = "dense"
      else
        matrix_type = "periodic_banded"
      end if
    else
      matrix_type = "banded"
    end if

    call sll_s_spline_matrix_new( self%matrix, matrix_type, self%n, kl, ku )
    call build_system( self, self%matrix )
    call self % matrix % factorize()

  end subroutine s_spline_1d_uniform__init

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for any kind of boundary conditions at xmin and xmax
  !> @param[in]    self   spline interpolation object
  !> @param[inout] matrix generic 'spline' matrix (dense/banded/periodic-banded)
  !-----------------------------------------------------------------------------
  subroutine build_system( self, matrix )
    class(sll_t_spline_1d_uniform), intent(in   ) :: self
    class(sll_c_spline_matrix)    , intent(inout) :: matrix

    integer  :: i,j,d,s
    integer  :: j0,d0
    integer  :: icell
    integer  :: order
    real(wp) :: x
    real(wp) :: x_offset
    real(wp) :: values(self%deg+1)
    real(wp), allocatable :: derivs(:,:)

    if ( any( [self%bc_xmin,self%bc_xmax] == sll_p_hermite ) ) then
      allocate ( derivs (0:self%deg/2, 1:self%deg+1) )
    end if

    ! Hermite boundary conditions at xmin, if any
    ! NOTE: When using uniform B-splines, the cell size is normalized between
    !       0 and 1, hence the i-th derivative is multiplied by dx^i.
    !       This is beneficial to the linear system, acting as a natural
    !       preconditioner by avoiding very large matrix entries, and hence
    !       reducing the amplification of round-off errors.
    if ( self%bc_xmin == sll_p_hermite ) then
      x        = self%xmin
      icell    = 1
      x_offset = 0.0_wp
      call sll_s_uniform_bsplines_eval_basis_and_n_derivs( self%deg, x_offset, self%nbc_xmin, derivs )
      do i = 1, self%nbc_xmin
        order = self%nbc_xmin-i+self%mod
        do j = 1, self%deg  ! iterate only to deg as last bspline is 0
          call matrix % set_element( i, j, derivs(order,j) )
        end do
      end do
      if ( self%even ) then
        call sll_s_uniform_bsplines_eval_basis( self%deg, x_offset, values )
        i = self%nbc_xmin
        do j = 1, self%deg  ! iterate only to deg as last bspline is 0
          call matrix % set_element( i, j, values(j) )
        end do
      end if
    end if

    ! Interpolation points
    ! TODO: if x_offset close to 0, set x_offset=0
    ! TODO: if x_offset close to 1, set x_offset=0 and icell=icell-1
    do i = self%nbc_xmin+1, self%n-self%nbc_xmax
      x = self%tau(i-self%nbc_xmin)
      call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, x_offset )
      call sll_s_uniform_bsplines_eval_basis( self%deg, x_offset, values )
      do s = 1, self%deg+1
        j = modulo(icell-self%offset-2+s,self%n)+1
        call matrix % set_element( i, j, values(s) )
      end do
    end do

    ! Hermite boundary conditions at xmax, if any
    ! NOTE: When using uniform B-splines, the cell size is normalized between
    !       0 and 1, hence the i-th derivative is multiplied by dx^i.
    !       This is beneficial to the linear system, acting as a natural
    !       preconditioner by avoiding very large matrix entries, and hence
    !       reducing the amplification of round-off errors.
    if ( self%bc_xmax == sll_p_hermite ) then
      x        = self%xmax
      icell    = self%ncells
      x_offset = 1.0_wp
      call sll_s_uniform_bsplines_eval_basis_and_n_derivs( self%deg, x_offset, self%nbc_xmax, derivs )
      do i = self%n-self%nbc_xmax+1, self%n
        order = i-(self%n-self%nbc_xmax+1)+self%mod
        j0 = self%n-self%deg
        d0 = 1
        do s = 1, self%deg
          j = j0 + s
          d = d0 + s
          call matrix % set_element( i, j, derivs(order,d) )
        end do
      end do
      if ( self%even ) then
        call sll_s_uniform_bsplines_eval_basis( self%deg, x_offset, values )
        i  = self%n-self%nbc_xmax+1
        j0 = self%n-self%deg
        d0 = 1
        do s = 1, self%deg
          j = j0 + s
          d = d0 + s
          call matrix % set_element( i, j, values(d) )
        end do
      end if
    end if

    if ( allocated( derivs ) ) deallocate( derivs )

  end subroutine build_system

!  !-----------------------------------------------------------------------------
!  !> @brief        Private subroutine for assembling and factorizing linear
!  !>               system needed for periodic spline interpolation
!  !> @param[inout] self spline 1D object
!  !-----------------------------------------------------------------------------
!  subroutine s_build_system_periodic( self )
!
!    class(sll_t_spline_1d_uniform), intent(inout) :: self
!
!    integer  :: j
!    integer  :: k
!    real(wp) :: offset
!    real(wp) :: values(self%deg+1)
!
!    ! Evaluate B-splines at interpolation points, noting that:
!    !   a) points are regularly spaced
!    !   b) basis is translationally invariant
!    ! For even degree interpolation points are cell midpoints
!    ! For odd  degree interpolation points are grid points
!    if (modulo(self%deg,2) == 0) then
!      offset = 0.5_wp
!    else
!      offset = 0.0_wp
!    end if
!    call sll_s_uniform_bsplines_eval_basis( self%deg, offset, values )
!
!    ! Allocate matrix Q for linear system:
!    !   - M is circulant banded matrix (sparse)
!    !   - Q is compressed form of A (circulant diagonals of M are rows of Q)
!    k = self%deg/2
!    allocate( self%q(2*k+1,self%n) )
!
!    ! Assemble matrix Q:
!    !  - M has constant elements on each diagonal
!    !  - Q has constant elements on each row
!    do j = 1, self%n
!      self%q(:,j) = values(1:2*k+1)
!    end do
!
!    ! Prepare solution of linear system by performing:
!    ! a) Block decomposition of M = [[A,B],[C,D]] with
!    !     - A (n-k)x(n-k) banded Toeplitz matrix (square, large)
!    !     - B (n-k)x(  k) low-rank Toeplitz (rectangular, empty in the center)
!    !     - C (  k)x(n-k) low-rank Toeplitz (rectangular, empty in the center)
!    !     - D (  k)x(  k) fully dense Toeplitz matrix (square, small)
!    ! b) LU factorization of matrix A
!    ! c) Calculation of Schur complement of A: H = D - C A^{-1} B
!    ! d) LU factorization of Schur complement H
!    call schur_complement_fac( self%schur, self%n, k, self%q )
!
!  end subroutine s_build_system_periodic
!
!  !-----------------------------------------------------------------------------
!  !> @brief        Private subroutine for assembling and factorizing linear
!  !>               system for spline interpolation with Hermite BCs
!  !> @param[inout] self spline 1D object
!  !-----------------------------------------------------------------------------
!  subroutine s_build_system_hermite( self )
!
!    class(sll_t_spline_1d_uniform), intent(inout) :: self
!
!    integer  :: nbc
!    integer  :: iflag
!    integer  :: i
!    integer  :: j
!    integer  :: ii
!    integer  :: jj
!    integer  :: offset
!    integer  :: k
!    integer  :: icell
!    real(wp) :: x
!    real(wp) :: x_offset
!
!    real(wp) :: values(self%deg+1)
!    real(wp) :: derivs(0:self%deg/2,1:self%deg+1)
!
!
!    ! number of boundary conditions needed depending on spline degree
!    k   = self%deg-1
!    nbc = self%deg/2
!
!    ! The matrix q is a banded matrix using the storage required by DGBTRF (LAPACK)
!    ! It has 2*k bands above diagonal, k bands below and the diagonal itself
!    ! k additional bands above diagonal needed for pivoting
!    ! The term A(ii,jj) of the full matrix is stored in q(ii-jj+2*k+1,jj)
!    ! The Bspline interpolation matrix at Greville points x_ii is
!    ! A(ii,jj) = B_jj(x_ii)
!    allocate( self%q(3*k+1,self%n) )
!    self%q = 0.0_wp
!
!    ! For even degree splines interpolation points are at cell midpoints
!    ! value of function on the boundary needed as additional boundary conditions
!    ! for odd degree splines  baundary is in the interpolation points
!    ! only derivative values are needed as boundary conditions
!    if (modulo( self%deg, 2 ) == 0) then ! spline degree even
!       offset = 1
!    else
!       offset = 0
!    end if
!
!    ! boundary conditions at xmin
!    ii       = 0
!    x        = self%xmin
!    icell    = 1
!    x_offset = 0.0_wp
!    call sll_s_uniform_bsplines_eval_basis_and_n_derivs( self%deg, x_offset, nbc, derivs )
!    ! When using uniform B-splines, the cell size is normalized between 0 and 1,
!    ! hence the i-th derivative must be divided by dx^i to scale back the result
!    do i = 1, nbc-offset
!      derivs(i,:) = derivs(i,:) / self%dx**i
!    end do
!    do i = 1, nbc
!       ! iterate only to k+1 as last bspline is 0
!       ii = ii + 1
!       do j = 1, k+1
!          jj = icell+j-1
!          self%q(ii-jj+2*k+1,jj) = derivs(i-offset,j)
!       end do
!    end do
!
!    ! interpolation points
!    do i = 1, self%n-2*nbc
!       ii = ii + 1
!       x = self%tau(i)
!       call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, x_offset )
!       call sll_s_uniform_bsplines_eval_basis( self%deg, x_offset, values )
!       do j=1,k+2
!          jj = icell+j-1
!          self%q(ii-jj+2*k+1,jj) = values(j)
!       end do
!    end do
!
!    ! boundary conditions at xmax
!    x        = self%xmax
!    icell    = self%n-self%deg
!    x_offset = 1.0_wp
!    call sll_s_uniform_bsplines_eval_basis_and_n_derivs( self%deg, x_offset, nbc, derivs )
!    ! When using uniform B-splines, the cell size is normalized between 0 and 1,
!    ! hence the i-th derivative must be divided by dx^i to scale back the result
!    do i = 1, nbc-offset
!      derivs(i,:) = derivs(i,:) / self%dx**i
!    end do
!    do i = 1, nbc
!       ii = ii + 1
!       do j = 2, k+2
!          jj = icell+j-1
!          self%q(ii-jj+2*k+1,jj) = derivs(i-offset,j)
!       end do
!    end do
!
!    ! Perform LU decomposition of matrix q with Lapack
!    allocate( self%ipiv (self%n) )
!    call dgbtrf( self%n, self%n, k, k, self%q, 3*k+1, self%ipiv, iflag )
!
!  end subroutine s_build_system_hermite
!
!  !-----------------------------------------------------------------------------
!  !> @brief        Private subroutine for assembling and factorizing linear
!  !>               system for spline interpolation at Greville points
!  !> @param[inout] self bspline object
!  !-----------------------------------------------------------------------------
!  subroutine s_build_system_greville( self )
!
!    class(sll_t_spline_1d_uniform), intent(inout) :: self
!
!    integer  :: iflag
!    integer  :: j
!    integer  :: k
!    integer  :: ii
!    integer  :: jj
!    integer  :: icell
!    real(wp) :: x
!    real(wp) :: offset
!    real(wp) :: values(self%deg+1)
!
!    ! The matrix q is a banded matrix using the storage required by banfac (De Boor)
!    ! It has k bands above diagonal, k bands below and the diagonal itself
!    ! The term A(ii,jj) of the full matrix is stored in q(ii-jj+k+1,jj)
!    ! The Bspline interpolation matrix at Greville points x_ii is
!    ! A(ii,jj) = B_jj(x_ii)
!    k = self%deg - 1
!    allocate( self%q(2*k+1,self%n) )
!    self%q = 0.0_f64
!
!    ! Treat i=1 separately
!    x     =  self%xmin
!    icell = 1
!    ii    = 1   ! first Greville point
!    call sll_s_uniform_bsplines_eval_basis( self%deg, 0.0_wp, values )
!    ! iterate only to k+1 as last bspline is 0
!    do j=1,k+1
!       jj = icell+j-1
!       self%q(ii-jj+k+1,jj) = values(j)
!    end do
!
!    do ii=2,self%n - 1
!       x = self%tau(ii)
!       call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, offset )
!       call sll_s_uniform_bsplines_eval_basis( self%deg, offset, values )
!       do j=1,k+2
!          jj = icell+j-1
!          self%q(ii-jj+k+1,jj) = values(j)
!       end do
!    end do
!
!    ! Treat i=self%n separately
!    x     = self%xmax
!    icell = self%n - self%deg
!    ii    = self%n  ! last Greville point
!    call sll_s_uniform_bsplines_eval_basis( self%deg, 1.0_wp, values )
!    ! iterate only to k as first bspline is 0
!    do j=2,k+2
!       jj = icell+j-1
!       self%q(ii-jj+k+1,jj) = values(j)
!    end do
!
!    ! Perform LU decomposition of matrix q
!    call banfac ( self%q, 2*k+1, self%n, k, k, iflag )
!
!  end subroutine s_build_system_greville

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

    character(len=*), parameter :: &
      this_sub_name = "spline_1d_non_uniform % compute_interpolant"

    integer :: i

    ! Special case: linear spline
    if (self%deg == 1) then
      self%bcoef(1:self%n) = gtau(1:self%n)
      ! Periodic only: "wrap around" coefficients onto extended array
      if (self%bc_xmin == sll_p_periodic) then
        self%bcoef(0)        = self%bcoef(self%n)
        self%bcoef(self%n+1) = self%bcoef(1)
      end if
      return
    end if

    ! Hermite boundary conditions at xmin, if any
    ! NOTE: When using uniform B-splines, the cell size is normalized between
    !       0 and 1, hence for consistency with the linear system, the i-th
    !       derivative provided by the user must be multiplied by dx^i
    if ( self%bc_xmin == sll_p_hermite ) then
      if ( present( derivs_xmin ) ) then
        self%bcoef(1:self%nbc_xmin) = &
                [(derivs_xmin(i)*self%dx**(i+self%mod-1), i=self%nbc_xmin,1,-1)]
      else
        SLL_ERROR(this_sub_name,"Hermite BC at xmin requires derivatives")
      end if
    end if

    ! interpolation points
    self%bcoef(self%nbc_xmin+1:self%n-self%nbc_xmax) = gtau(:)

    ! Hermite boundary conditions at xmax, if any
    ! NOTE: When using uniform B-splines, the cell size is normalized between
    !       0 and 1, hence for consistency with the linear system, the i-th
    !       derivative provided by the user must be multiplied by dx^i
    if ( self%bc_xmax == sll_p_hermite ) then
      if ( present( derivs_xmax ) ) then
        self%bcoef(self%n-self%nbc_xmax+1:self%n) = &
                   [(derivs_xmax(i)*self%dx**(i+self%mod-1), i=1,self%nbc_xmax)]
      else
        SLL_ERROR(this_sub_name,"Hermite BC at xmax requires derivatives")
      end if
    end if

    ! Solve linear system and compute coefficients
    call self % matrix % solve_inplace( self%bcoef(1:self%n) )

    ! Periodic only: "wrap around" coefficients onto extended array
    if (self%bc_xmin == sll_p_periodic) then
      associate( n => self%n, g => 1+self%deg/2 )
        self%bcoef(1-g:0)   = self%bcoef(n-g+1:n)
        self%bcoef(n+1:n+g) = self%bcoef(1:g)
      end associate
    end if

  end subroutine s_spline_1d_uniform__compute_interpolant

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d_uniform__eval( self, x ) result( y )

    class(sll_t_spline_1d_uniform), intent(in) :: self
    real(wp)                      , intent(in) :: x
    real(wp) :: y

    integer  :: icell
    real(wp) :: offset
    real(wp) :: values(self%deg+1)

    call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, offset )
    call sll_s_uniform_bsplines_eval_basis( self%deg, offset, values )

    associate( a => icell-self%offset, &
               b => icell-self%offset+self%deg )
      y = dot_product( self%bcoef(a:b), values )
    end associate

  end function f_spline_1d_uniform__eval

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d_uniform__eval_deriv( self, x ) result( y )

    class(sll_t_spline_1d_uniform), intent(in) :: self
    real(wp)                      , intent(in) :: x
    real(wp) :: y

    integer  :: icell
    real(wp) :: offset
    real(wp) :: values(self%deg+1)

    call s_get_icell_and_offset( self%xmin, self%xmax, self%ncells, x, icell, offset )
    call sll_s_uniform_bsplines_eval_deriv( self%deg, offset, values )

    associate( a => icell-self%offset, &
               b => icell-self%offset+self%deg )
      y = dot_product( self%bcoef(a:b), values )
    end associate

    ! When using uniform B-splines, the cell size is normalized between 0 and 1,
    ! hence the derivative must be divided by dx to scale back the result.
    y = y * self%inv_dx

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

    ! Free matrix and linear solver
    if (allocated( self%matrix )) then
      call self % matrix % free()
      deallocate( self % matrix )
    end if

  end subroutine s_spline_1d_uniform__free

end module sll_m_spline_1d_uniform
