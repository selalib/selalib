!> @ingroup splines
!> @brief Arbitrary degree spline interpolation on a non-uniform grid
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_1d_non_uniform
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite , &
    sll_p_greville, &
    sll_p_open
!    sll_p_open    , &
!    sll_p_mirror

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

  use sll_m_spline_matrix, only: &
    sll_c_spline_matrix, &
    sll_s_spline_matrix_new

  implicit none

  public :: sll_t_spline_1d_non_uniform

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Allowed boundary conditions
  integer, parameter :: allowed_bcs(*) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  !> 1D spline interpolation on non-uniform grid
  type, extends(sll_c_spline_1d) :: sll_t_spline_1d_non_uniform

    integer               :: deg      ! spline degree (= order of piecewise polynomial)
    integer               :: n        ! dimension of spline space
    integer               :: bc_xmin  ! boundary condition type at x=xmin
    integer               :: bc_xmax  ! boundary condition type at x=xmax
    integer               :: offset   ! needed for periodic spline evaluation
    type (sll_t_bsplines) :: bsp      ! basis functions (B-splines)

    real(wp), private :: xmin     ! left  boundary coordinate
    real(wp), private :: xmax     ! right boundary coordinate
    integer , private :: nbc_xmin ! number of boundary conditions (derivatives) at x=xmin
    integer , private :: nbc_xmax ! number of boundary conditions (derivatives) at x=xmax
    integer , private :: ncells   ! number of cells
    integer , private :: mod      ! result of modulo(deg,2): 0 if deg even, 1 if deg odd

    real(wp), private :: dx       ! cell size (averaged over non-uniform domain)

    real(wp), allocatable :: bcoef(:) ! B-splines' coefficients
    real(wp), allocatable :: tau(:)   ! interpolation points

    ! Polymorphic matrix object to store and solve linear system for interpolation
    class(sll_c_spline_matrix), allocatable, private :: matrix

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
    integer :: i,kl,ku
    integer :: basis_bc_xmin
    integer :: basis_bc_xmax

    real(wp), allocatable :: temp_knots(:)

    character(len=*), parameter :: this_sub_name = "spline_1d_non_uniform % init"
    character(len=32) :: matrix_type

    ! Sanity checks
    SLL_ASSERT( degree >= 1 )
    SLL_ASSERT( size( breaks ) >= 2 )
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )

    if ( any( [bc_xmax,bc_xmin] == sll_p_periodic ) .and. bc_xmin /= bc_xmax) then
      SLL_ERROR(this_sub_name,"Periodic BCs cannot be mixed with Hermite or Greville BCs")
    end if

    self%bc_xmin = bc_xmin
    self%bc_xmax = bc_xmax

    self%ncells = size(breaks) - 1

    self%xmin = breaks(1)
    self%xmax = breaks(self%ncells+1)
    self%dx   = (self%xmax-self%xmin)/self%ncells

    ! set first attributes
    if ( self%bc_xmin == sll_p_periodic) then
      self%n = self%ncells    ! dimension of periodic spline space
      self%offset = degree/2  ! offset needed for periodic spline evaluation
    else
      self%n = self%ncells + degree ! dimension of non periodic spline space
      self%offset = 0
    end if

    self%deg = degree
    self%mod = modulo( degree, 2 )

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

    ! Allocate array of spline coefficients
    ! in case of periodic BCs, a larger array of coefficients is used in order
    ! to avoid a loop with calls to the "mod( , )" function at evaluation.
    associate( n => self%n, g => 1+self%deg/2 )
      if (self%bc_xmin == sll_p_periodic) then
        allocate( self%bcoef(1-g:n+g) )
      else
        allocate( self%bcoef(1:n) )
      end if
    end associate

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
    else if ( self%mod == 1 ) then
      self%tau(1)    = self%xmin
      self%tau(ntau) = self%xmax
    end if

    ! Special case: linear spline (no need for matrix assembly)
    if (self%deg == 1) return

    !---------------------------------------------------------------------------
    ! Assemble dense matrix (B_j(tau(i))) for spline interpolation
    !---------------------------------------------------------------------------

    ! FIXME: In Hermite case ku and kl computed in general case when derivatives
    !        of B-splines do not vanish at boundary
    select case( self%bc_xmin )
      case ( sll_p_periodic ); ku = ( self%deg + 1 ) / 2
      case ( sll_p_hermite  ); ku = max( (self%deg+1)/2, self%deg-1 )
      case ( sll_p_greville ); ku = self%deg
    end select

    select case( self%bc_xmax )
      case ( sll_p_periodic ); kl = ( self%deg + 1 ) / 2
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

  end subroutine s_spline_1d_non_uniform__init

  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for any kind of boundary conditions at xmin and xmax
  !> @param[in]    self   spline interpolation object
  !> @param[inout] matrix generic 'spline' matrix (dense/banded/periodic-banded)
  !-----------------------------------------------------------------------------
  subroutine build_system( self, matrix )
    class(sll_t_spline_1d_non_uniform), intent(in   ) :: self
    class(sll_c_spline_matrix)        , intent(inout) :: matrix

    integer  :: i,j,d,s
    integer  :: j0,d0
    integer  :: icell
    integer  :: order
    real(wp) :: x
    real(wp) :: values(self%deg+1)
    real(wp), allocatable :: derivs(:,:)

    if ( any( [self%bc_xmin,self%bc_xmax] == sll_p_hermite ) ) then
      allocate ( derivs (0:self%deg/2, 1:self%deg+1) )
    end if

    ! Hermite boundary conditions at xmin, if any
    if ( self%bc_xmin == sll_p_hermite ) then
      x = self%xmin
      icell = 1
      call sll_s_bsplines_eval_basis_and_n_derivs( self%bsp, icell, x , self%nbc_xmin, derivs )

      ! In order to improve the condition number of the matrix, we normalize all
      ! derivatives by multiplying the i-th derivative by dx^i
      associate( h => [(self%dx**i, i=1, ubound(derivs,1))] )
        do j = lbound(derivs,2), ubound(derivs,2)
          derivs(1:,j) = derivs(1:,j) * h(1:)
        end do
      end associate

      do i = 1, self%nbc_xmin
        ! iterate only to deg as last bspline is 0
        order = self%nbc_xmin-i+self%mod
        do j = 1, self%deg
          call matrix % set_element( i, j, derivs(order,j) )
        end do
      end do

    end if

    ! Interpolation points
    do i = self%nbc_xmin+1, self%n-self%nbc_xmax
      x = self%tau(i-self%nbc_xmin)
      icell = sll_f_find_cell( self%bsp, x )
      call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )
      do s = 1, self%deg+1
        j = modulo(icell-self%offset-2+s,self%n)+1
        call matrix % set_element( i, j, values(s) )
      end do
    end do

    ! Hermite boundary conditions at xmax, if any
    if ( self%bc_xmax == sll_p_hermite ) then
      x = self%xmax
      icell = self%ncells
      call sll_s_bsplines_eval_basis_and_n_derivs( self%bsp, icell, x, self%nbc_xmax, derivs )

      ! In order to improve the condition number of the matrix, we normalize all
      ! derivatives by multiplying the i-th derivative by dx^i
      associate( h => [(self%dx**i, i=1, ubound(derivs,1))] )
        do j = lbound(derivs,2), ubound(derivs,2)
          derivs(1:,j) = derivs(1:,j) * h(1:)
        end do
      end associate

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

    end if

    if ( allocated( derivs ) ) deallocate( derivs )

  end subroutine build_system

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
    real(wp)                          , intent(in   ) :: gtau(:)
    real(wp),                 optional, intent(in   ) :: derivs_xmin(:)
    real(wp),                 optional, intent(in   ) :: derivs_xmax(:)

    character(len=*), parameter :: &
      this_sub_name = "spline_1d_non_uniform % compute_interpolant"

    integer :: i

    SLL_ASSERT( size(gtau) == self%n-self%nbc_xmin-self%nbc_xmax )

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
    ! NOTE: For consistency with the linear system, the i-th derivative
    !       provided by the user must be multiplied by dx^i
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
    ! NOTE: For consistency with the linear system, the i-th derivative
    !       provided by the user must be multiplied by dx^i
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
    real(wp) :: values(self%deg+1)

    icell = sll_f_find_cell( self%bsp, x )
    call sll_s_bsplines_eval_basis( self%bsp, icell, x, values )

    associate( a => icell-self%offset, &
               b => icell-self%offset+self%deg )
      y = dot_product( self%bcoef(a:b), values )
    end associate

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
    real(wp) :: values(self%deg+1)

    icell = sll_f_find_cell( self%bsp, x )
    call sll_s_bsplines_eval_deriv( self%bsp, icell, x, values )

    associate( a => icell-self%offset, &
               b => icell-self%offset+self%deg )
      y = dot_product( self%bcoef(a:b), values )
    end associate

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

    ! Deallocate arrays
    deallocate( self%bcoef )
    deallocate( self%tau   )

    ! Free B-splines' memory
    call sll_s_bsplines_free( self%bsp )

    ! Free matrix and linear solver
    if (allocated( self % matrix )) then
      call self % matrix % free()
      deallocate( self % matrix )
    end if

  end subroutine s_spline_1d_non_uniform__free

end module sll_m_spline_1d_non_uniform
