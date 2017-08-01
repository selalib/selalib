!> @ingroup splines
!> Implements arbitrary degree bspline interpolation on a uniform grid
!> given a B-Spline object from sll_m_bsplines

module sll_m_spline_1d_non_uniform
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

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

  public :: sll_t_spline_1d_non_uniform

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Allowed boundary conditions
  integer, parameter :: allowed_bcs(*) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  ! TODO: use this to employ banded matrix q instead of dense matrix a
  integer, parameter :: m = 0 ! additional diagonals above and below

  !> @brief
  !> basic type for one-dimensional B-spline interpolation.
  type, extends(sll_c_spline_1d) :: sll_t_spline_1d_non_uniform

    type (sll_t_bsplines) :: bsp
    integer               :: n        ! dimension of spline space
    integer               :: deg      ! degree of spline
    real(wp), allocatable :: tau(:)   ! interpolation points
    real(wp), allocatable :: bcoef(:) ! bspline coefficients
    integer               :: bc_type  ! boundary condition
    real(wp), allocatable :: q(:,:)   ! triangular factorization
    real(wp), allocatable :: bsdx(:,:)
    integer , allocatable :: ipiv(:)
    integer               :: offset
    type(schur_complement_solver) :: schur

    real(wp) :: xmin
    real(wp) :: xmax
    integer  :: bc_xmin  ! boundary condition type at x=xmin
    integer  :: bc_xmax  ! boundary condition type at x=xmax
    integer  :: nbc_xmin ! number of boundary conditions at x=xmin
    integer  :: nbc_xmax ! number of boundary conditions at x=xmax
    integer  :: ncells   ! number of cells
    logical  :: even     ! true if deg even, false if deg odd
    integer  :: mod      ! result of modulo(deg,2): 0 if deg even, 1 if deg odd
    real(wp), allocatable :: a(:,:) ! dense matrix (alternative to q)

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

  !-----------------------------------------------------------------------------
  function f_spline_1d_non_uniform__get_coeff( self ) result( ptr )

    class(sll_t_spline_1d_non_uniform), target, intent(in) :: self
    real(wp), pointer :: ptr(:)

    ptr => self%bcoef

  end function f_spline_1d_non_uniform__get_coeff

  !-----------------------------------------------------------------------------
  subroutine s_spline_1d_non_uniform__get_interp_points( self, tau )

    class(sll_t_spline_1d_non_uniform), intent(in   ) :: self
    real(wp),              allocatable, intent(  out) :: tau(:)

    SLL_ASSERT( allocated( self%tau ) )
    allocate( tau(size(self%tau)), source=self%tau )

  end subroutine s_spline_1d_non_uniform__get_interp_points

  !-----------------------------------------------------------------------------
  !> @brief     Constructor for sll_t_bspline_1d object
  !> @param[in] degree   Spline degree
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
    integer                           , intent(in   ) :: degree
    real(wp)                          , intent(in   ) :: breaks(:)
    integer                           , intent(in   ) :: bc_xmin
    integer                           , intent(in   ) :: bc_xmax

    integer :: ntau
    integer :: i
    integer :: bc_type
    integer :: basis_bc_xmin
    integer :: basis_bc_xmax

    real(wp), allocatable :: temp_knots(:)

    character(len=*), parameter :: this_sub_name = "spline_1d_non_uniform % init"

    ! Sanity checks
    SLL_ASSERT( degree >= 1 )
    SLL_ASSERT( size( breaks ) >= 2 )
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )

    if ( any( [bc_xmax,bc_xmin] == sll_p_periodic ) .and. bc_xmin /= bc_xmax) then
      SLL_ERROR(this_sub_name,"Periodic BCs cannot be mixed with Hermite or Greville BCs")
    end if

    bc_type = bc_xmin

    self%bc_xmin = bc_xmin
    self%bc_xmax = bc_xmax

    self%ncells = size(breaks) - 1

    self%xmin = breaks(1)
    self%xmax = breaks(self%ncells+1)

    ! set first attributes
    self%bc_type = bc_type
    if (bc_type == sll_p_periodic) then
      self%n = self%ncells    ! dimension of periodic spline space
      self%offset = degree/2  ! offset needed for periodic spline evaluation
    else
      self%n = self%ncells + degree ! dimension of non periodic spline space
      self%offset = 0
    end if

    self%deg = degree

    self%mod = modulo( degree, 2 )

    if ( self%mod == 0 ) then
      self%even = .true.
    else
      self%even = .false.
    end if

    self%nbc_xmin = merge( degree/2, 0, self%bc_xmin == sll_p_hermite )
    self%nbc_xmax = merge( degree/2, 0, self%bc_xmax == sll_p_hermite )

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
    allocate( self%bsdx (degree/2+1,degree+1) )
    allocate( self%ipiv (self%n) )

    !---------------------------------------------------------------------------
    ! Determine array tau of interpolation points
    !---------------------------------------------------------------------------

    ! Determine size of tau and allocate tau
    if ( self%bc_xmin == sll_p_periodic ) then
      ntau = self%ncells
    else
      ntau = self%ncells + degree - self%nbc_xmin - self%nbc_xmax
    end if
    allocate( self%tau( 1 : ntau ) )

    ! Array of temporary knots needed to compute interpolation points
    ! using Greville-style averaging: tau(i) = average(temp_knots(i+1-degree:i))
    allocate( temp_knots( 2-degree : ntau ) )

    if ( self%bc_xmin == sll_p_periodic ) then

      associate( k => degree/2 )
        temp_knots(:) = self%bsp%knots(2-degree+k:ntau+k)
      end associate

    else

      associate( r => 2-degree, s => -self%nbc_xmin )
        select case (bc_xmin)
        case (sll_p_greville); temp_knots(r:s) = breaks(1)
        case (sll_p_hermite ); temp_knots(r:s) = 2.0_wp*breaks(1) - breaks(2+s-r:2:-1)
        end select
      end associate

      associate( r => -self%nbc_xmin+1, s => -self%nbc_xmin+1+self%ncells )
        temp_knots(r:s) = breaks(:)
      end associate

      associate( r => -self%nbc_xmin+1+self%ncells+1, s => ntau, n => self%ncells )
        select case (bc_xmax)
        case (sll_p_greville); temp_knots(r:s) = breaks(n+1)
        case (sll_p_hermite ); temp_knots(r:s) = 2.0_wp*breaks(n+1) - breaks(n:n+r-s:-1)
        end select
      end associate

    end if

    ! Compute interpolation points using Greville-style averaging
    associate( inv_deg => 1.0_wp / real( degree, wp ) )
      do i = 1, ntau
        self%tau(i) = sum( temp_knots(i+1-degree:i) ) * inv_deg
      end do
    end associate

    ! Periodic case: apply periodic BCs to interpolation points
    if ( self%bc_xmin == sll_p_periodic ) then
      self%tau(:) = modulo( self%tau(:)-self%xmin, self%xmax-self%xmin ) + self%xmin

    ! Non-periodic case, odd degree: fix round-off issues
    else if ( .not. self%even ) then
      self%tau(1)    = self%xmin
      self%tau(ntau) = self%xmax
    end if

    ! Special case: linear spline (no need for matrix assembly)
    if (self%deg == 1) then
      allocate( self%q(0,0) )
      return
    end if

    ! Assemble dense matrix (B_j(tau(i))) for spline interpolation
    call build_system_dense( self )

  end subroutine s_spline_1d_non_uniform__init

  subroutine build_system_dense( self )
    class(sll_t_spline_1d_non_uniform), intent(inout) :: self

    integer  :: i,j,s
    integer  :: icell
    integer  :: iflag
    integer  :: order
    real(wp) :: x
    real(wp) :: values(self%deg+1)
    real(wp), allocatable :: derivs(:,:)

!    ! DEBUG
!    real(wp), allocatable :: at(:,:)
!    character(len=20) :: fmt

    allocate( self%a(self%n,self%n) )
    self%a(:,:) = 0.0_wp

!    ! DEBUG
!    allocate( at(self%n,self%n) )

    if ( any( [self%bc_xmin,self%bc_xmax] == sll_p_hermite ) ) then
      allocate ( derivs (0:self%deg/2, 1:self%deg+1) )
    end if

    ! Hermite boundary conditions at xmin, if any
    if ( self%bc_xmin == sll_p_hermite ) then
      x = self%xmin
      icell = 1
      call sll_s_bsplines_eval_basis_and_n_derivs( self%bsp, icell, x , self%nbc_xmin, derivs )
      do i = 1, self%nbc_xmin
        ! iterate only to deg as last bspline is 0
        order = self%nbc_xmin-i+self%mod
        self%a(i,1:self%deg) = derivs(order,1:self%deg)
      end do
      if ( self%even ) then
        call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
        self%a(self%nbc_xmin,1:self%deg) = values(1:self%deg)
      end if
    end if

    ! interpolation points
    do i = self%nbc_xmin+1, self%n-self%nbc_xmax
      x = self%tau(i-self%nbc_xmin)
      icell = sll_f_find_cell( self%bsp, x )
      call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
      do s = 1, self%deg+1
        j = modulo(icell-self%offset-2+s,self%n)+1
        self%a(i,j) = values(s)
      end do
    end do

    ! Hermite boundary conditions at xmax, if any
    if ( self%bc_xmax == sll_p_hermite ) then
      x = self%xmax
      icell = self%ncells
      call sll_s_bsplines_eval_basis_and_n_derivs( self%bsp, icell, x , self%nbc_xmax, derivs )
      do i = self%n-self%nbc_xmax+1, self%n
        order = i-(self%n-self%nbc_xmax+1)+self%mod
        self%a(i,self%n-self%deg+1:self%n) = derivs(order,2:self%deg+1)
      end do
      if ( self%even ) then
        call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
        self%a(self%n-self%nbc_xmax+1,self%n-self%deg+1:self%n) = values(2:self%deg+1)
      end if
    end if

!    ! DEBUG
!    do j = 1, self%n
!      do i = 1, self%n
!        at(i,j) = self%a(j,i)
!      end do
!    end do
!
!    write(*,*)
!    write(*,'(*(f8.4))') self%tau
!    write(*,*)
!
!    write(fmt,"('(',i0,'es9.1)')") self%n
!    write(*,*)
!    write(*,fmt) at
!    write(*,*)

    call dgetrf( self%n, self%n, self%a, self%n, self%ipiv, iflag )

    if ( allocated( derivs ) ) deallocate( derivs )

  end subroutine build_system_dense

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
    !allocate( self%q(2*k+1,self%n) )
    allocate( self%q(2*(k+m)+1,self%n) )

    self%q(:,:) = 0.0_wp

    ! evaluate bsplines at interpolation points
    ! and assemble matrix q for linear system
    do i = 1, self%n
      x = self%tau(i)
      icell = sll_f_find_cell( self%bsp, x )
      call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
      do s = 1, 2*k+1
        j = modulo(icell-k-2+s,self%n)+1
        !self%q(2*k+2-s,j) = values(s)
        self%q(m+2*k+2-s,j) = values(s)
      end do
    end do

    ! Perform LU decomposition of matrix q
    call schur_complement_fac( self%schur, self%n, k+m, self%q )

  end subroutine build_system_periodic

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for spline interpolation at Greville points
  !> @param[inout] self bspline object
  !-----------------------------------------------------------------------------
  subroutine build_system_greville( self )
    class(sll_t_spline_1d_non_uniform), intent(inout) :: self

    integer  :: j, k
    integer  :: ii, jj
    integer  :: icell
    integer  :: iflag
    real(wp) :: x
    real(wp) :: values(self%deg+1)

    ! The matrix q is a banded matrix using the storage required by banfac (De Boor)
    ! It has k bands above diagonal, k bands below and the diagonal itself
    ! The term A(ii,jj) of the full matrix is stored in q(ii-jj+k+1,jj)
    ! The Bspline interpolation matrix at Greville points x_ii is
    ! A(ii,jj) = B_jj(x_ii)
    k = self%deg - 1
    allocate( self%q(2*k+1,self%n) )
    self%q = 0.0_wp

    ! Treat i=1 separately
    x =  self%bsp%xmin
    icell = 1
    ii = 1   ! first Greville point
    call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
    ! iterate only to k+1 as last bspline is 0
    do j=1,k+1
      jj = icell+j-1
      self%q(ii-jj+k+1,jj) = values(j)
    end do

    do ii=2,self%n - 1
      x = self%tau(ii)
      icell = sll_f_find_cell(self%bsp, x )
      call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
      do j=1,k+2
        jj = icell+j-1
        self%q(ii-jj+k+1,jj) = values(j)
      end do
    end do

    ! Treat i=self%n separately
    x =  self%bsp%xmax
    icell = self%n - self%deg
    ii = self%n  ! last Greville point
    call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
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
    real(wp) :: values(self%deg+1)

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
    ! for odd degree splines boundary is in the interpolation points
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
!  subroutine s_spline_1d_non_uniform__compute_interpolant( self, gtau, derivs_xmin, derivs_xmax )
!    class(sll_t_spline_1d_non_uniform), intent(inout) :: self
!    real(wp)                          , intent(in   ) :: gtau(:)
!    real(wp),                 optional, intent(in   ) :: derivs_xmin(:)
!    real(wp),                 optional, intent(in   ) :: derivs_xmax(:)
!
!    integer :: ncond
!    integer :: k
!    integer :: iflag
!
!    ! Special case: linear spline
!    if (self%deg == 1) then
!      self%bcoef(:) = gtau(1:self%n)
!      return
!    end if
!
!    ! Degree > 1: boundary conditions matter
!    select case (self%bc_type)
!
!    case(sll_p_periodic)
!
!      self%bcoef = gtau(1:self%n)
!      call schur_complement_slv( self%schur, self%n, self%deg/2+m, self%q,  self%bcoef )
!
!    case (sll_p_greville)
!
!      self%bcoef = gtau
!      call banslv( self%q, 2*self%deg-1, self%n, self%deg-1, self%deg-1, self%bcoef )
!
!    case (sll_p_hermite)
!
!      ! number of needed conditions at boundary
!      ncond = self%deg/2
!
!      if (present(derivs_xmin)) then
!        self%bcoef(1:ncond) = derivs_xmin(ncond:1:-1)
!      else  ! set needed boundary values to 0
!        self%bcoef(1:ncond) = 0.0_wp
!      end if
!
!      self%bcoef(ncond+1:self%n-ncond) = gtau(1:self%n-2*ncond)
!
!      if (present(derivs_xmax)) then
!        self%bcoef(self%n-ncond+1:self%n) = derivs_xmax(1:ncond)
!      else ! set needed boundary values to 0
!        self%bcoef(self%n-ncond+1:self%n) = 0.0_wp
!      end if
!
!      k = self%deg-1
!      call dgbtrs( 'N', self%n, k, k, 1, self%q, 3*k+1, self%ipiv, self%bcoef, self%n, iflag )
!
!    end select
!
!  end subroutine s_spline_1d_non_uniform__compute_interpolant

  subroutine s_spline_1d_non_uniform__compute_interpolant( self, gtau, derivs_xmin, derivs_xmax )
    class(sll_t_spline_1d_non_uniform), intent(inout) :: self
    real(wp)                          , intent(in   ) :: gtau(:)
    real(wp),                 optional, intent(in   ) :: derivs_xmin(:)
    real(wp),                 optional, intent(in   ) :: derivs_xmax(:)

    integer :: iflag

    character(len=*), parameter :: &
      this_sub_name = "spline_1d_non_uniform % compute_interpolant"

    ! Special case: linear spline
    if (self%deg == 1) then
      self%bcoef(:) = gtau(1:self%n)
      return
    end if

    ! Hermite boundary conditions at xmin, if any
    if ( self%bc_xmin == sll_p_hermite ) then
      if ( present( derivs_xmin ) ) then
        self%bcoef(1:self%nbc_xmin) = derivs_xmin(self%nbc_xmin:1:-1)
      else
        SLL_ERROR(this_sub_name,"Hermite BC at xmin requires derivatives")
      end if
    end if

    ! interpolation points
    self%bcoef(self%nbc_xmin+1:self%n-self%nbc_xmax) = gtau(:)

    ! Hermite boundary conditions at xmax, if any
    if ( self%bc_xmax == sll_p_hermite ) then
      if ( present( derivs_xmax ) ) then
        self%bcoef(self%n-self%nbc_xmax+1:self%n) = derivs_xmax(1:self%nbc_xmax)
      else
        SLL_ERROR(this_sub_name,"Hermite BC at xmax requires derivatives")
      end if
    end if

    ! Solve linear system and compute coefficients
    call dgetrs( 'N', self%n, 1, self%a, self%n, self%ipiv, self%bcoef, self%n, iflag )

  end subroutine s_spline_1d_non_uniform__compute_interpolant

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate spline S at given point x
  !> @param[in]  self spline object
  !> @param[in]  x    evaluation point
  !> @return     y    spline value y=S(x)
  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_1d_non_uniform__eval( self, x ) result( y )
    class(sll_t_spline_1d_non_uniform), intent(in) :: self
    real(wp)                          , intent(in) :: x
    real(wp) :: y

    integer  :: icell
    integer  :: j, ib
    real(wp) :: values(self%deg+1)

    icell = sll_f_find_cell( self%bsp, x )
    call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
    y = 0.0_wp
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
    real(wp)                          , intent(in) :: x
    real(wp) :: y

    integer  :: icell
    integer  :: j, ib
    real(wp) :: values(self%deg+1)

    icell = sll_f_find_cell( self%bsp, x )
    call sll_s_bsplines_eval_deriv( self%bsp, icell, x, values )
    y = 0.0_wp
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
    real(wp)                          , intent(in   ) :: x(:)
    real(wp)                          , intent(  out) :: y(size(x))

    integer :: i

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
    real(wp)                          , intent(in   ) :: x(:)
    real(wp)                          , intent(  out) :: y(size(x))

    integer :: i

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
    deallocate( self%bcoef )
    deallocate( self%tau   )
    deallocate( self%bsdx  )
    if (allocated( self%q )) deallocate( self%q )

    ! free attribute objects
    call sll_s_bsplines_free( self%bsp )
    call schur_complement_free( self%schur )

  end subroutine s_spline_1d_non_uniform__free

end module sll_m_spline_1d_non_uniform
