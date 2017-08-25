!> @ingroup splines
!> Implements arbitrary degree bspline interpolation on a uniform grid
!> given a B-Spline object from sll_m_bsplines
!>
!> @author Eric Sonnendrücker - IPP Garching
!> @author Yaman Güçlü        - IPP Garching
!> @author Edoardo Zoni       - IPP Garching

module sll_m_spline_2d_uniform

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: &
    f64

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite, &
    sll_p_greville

  use sll_m_bsplines, only: &
    sll_s_uniform_bsplines_eval_basis, &
    sll_s_uniform_bsplines_eval_deriv

  use sll_m_spline_2d_base, only: &
    sll_c_spline_2d, &
    sll_t_spline_2d_boundary_data

  use sll_m_spline_1d_uniform, only: &
    sll_t_spline_1d_uniform

  implicit none

  public :: &
    sll_t_spline_2d_uniform

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> @brief
  !> basic type for two-dimensional B-spline data.
  !> @details
  !> treated as an opaque type. No access to its internals is directly allowed.
  type, extends(sll_c_spline_2d) :: sll_t_spline_2d_uniform

    integer ::     deg(2)
    integer ::       n(2)
    integer :: bc_xmin(2)
    integer :: bc_xmax(2)

    type(sll_t_spline_1d_uniform) :: bs1
    type(sll_t_spline_1d_uniform) :: bs2
    real(wp), allocatable         :: bcoef(:,:)
    real(wp)                      :: inv_dx(2)

    real(wp), allocatable, private :: bwork(:,:)
    integer              , private :: nbc_xmin(2)
    integer              , private :: nbc_xmax(2)

  contains
    ! Constructor
    procedure :: init                => s_spline_2d_uniform__init

    ! Abstract interface
    procedure :: free                => s_spline_2d_uniform__free
    procedure :: compute_interpolant => s_spline_2d_uniform__compute_interpolant
    procedure :: eval                => f_spline_2d_uniform__eval
    procedure :: eval_deriv_x1       => f_spline_2d_uniform__eval_deriv_x1
    procedure :: eval_deriv_x2       => f_spline_2d_uniform__eval_deriv_x2
    procedure :: eval_array          => s_spline_2d_uniform__eval_array
    procedure :: eval_array_deriv_x1 => s_spline_2d_uniform__eval_array_deriv_x1
    procedure :: eval_array_deriv_x2 => s_spline_2d_uniform__eval_array_deriv_x2
    procedure :: get_interp_points   => s_spline_2d_uniform__get_interp_points

  end type sll_t_spline_2d_uniform

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief Initialize a 2D spline interpolation object
  !-----------------------------------------------------------------------------
  subroutine s_spline_2d_uniform__init( &
    self   , &
    degree , &
    ncells , &
    xmin   , &
    xmax   , &
    bc_xmin, &
    bc_xmax )

    class(sll_t_spline_2d_uniform), intent(  out) :: self
    integer                       , intent(in   ) :: degree (2)
    integer                       , intent(in   ) :: ncells (2)
    real(wp)                      , intent(in   ) :: xmin   (2)
    real(wp)                      , intent(in   ) :: xmax   (2)
    integer                       , intent(in   ) :: bc_xmin(2)
    integer                       , intent(in   ) :: bc_xmax(2)

    call self%bs1%init( degree(1), ncells(1), xmin(1), xmax(1), bc_xmin(1), bc_xmax(1) )
    call self%bs2%init( degree(2), ncells(2), xmin(2), xmax(2), bc_xmin(2), bc_xmax(2) )

    associate( n1 => self%bs1%n, &
               n2 => self%bs2%n, &
               g1 => merge( 1+degree(1)/2, 0, bc_xmin(1) == sll_p_periodic ), &
               g2 => merge( 1+degree(2)/2, 0, bc_xmin(2) == sll_p_periodic ) )

      ! Save data into type
      self%deg     = degree
      self%n       = [n1,n2]
      self%bc_xmin = bc_xmin
      self%bc_xmax = bc_xmax
      self%inv_dx  = [real(ncells,wp)/(xmax-xmin)]

      ! Allocate array of spline coefficients, and work array for interpolation.
      ! in case of periodic BCs, a larger array of coefficients is used in order
      ! to avoid a loop with calls to the "mod( , )" function at evaluation.
      allocate( self%bcoef(1-g1:n1+g1,1-g2:n2+g2) );  self%bcoef = 0.0_f64
      allocate( self%bwork(1-g2:n2+g2,1-g1:n1+g1) );  self%bwork = 0.0_f64

      ! Calculate number of additional boundary data on each side of domain
      ! (i.e. derivatives for Hermite BC)
      self%nbc_xmin = merge( degree/2, 0, bc_xmin == sll_p_hermite )
      self%nbc_xmax = merge( degree/2, 0, bc_xmax == sll_p_hermite )

    end associate

  end subroutine s_spline_2d_uniform__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_2d_uniform__compute_interpolant( &
    self, &
    gtau, &
    boundary_data )

    class(sll_t_spline_2d_uniform)     , intent(inout)           :: self
    real(wp)                           , intent(in   )           :: gtau(:,:)
    type(sll_t_spline_2d_boundary_data), intent(in   ), optional :: boundary_data

    character(len=*), parameter :: &
      this_sub_name = "s_spline_2d_uniform__compute_interpolant"

    integer :: i1, i2

    ! User must provide derivatives at boundary in case of Hermite BC
    if (.not. present(boundary_data)) then
      if (self%bc_xmin(1) == sll_p_hermite) then
        SLL_ERROR( this_sub_name, "Hermite BC at x1_min requires derivatives" )
      else if (self%bc_xmax(1) == sll_p_hermite) then
        SLL_ERROR( this_sub_name, "Hermite BC at x1_max requires derivatives" )
      else if (self%bc_xmin(2) == sll_p_hermite) then
        SLL_ERROR( this_sub_name, "Hermite BC at x2_min requires derivatives" )
      else if (self%bc_xmax(2) == sll_p_hermite) then
        SLL_ERROR( this_sub_name, "Hermite BC at x2_max requires derivatives" )
      end if
    end if

    associate( w  => self%bcoef       , &
               wt => self%bwork       , &
               n1 => self%n(1)        , &
               a1 => self%nbc_xmin(1) , &
               b1 => self%nbc_xmax(1) , &
               n2 => self%n(2)        , &
               a2 => self%nbc_xmin(2) , &
               b2 => self%nbc_xmax(2) )

    SLL_ASSERT( all( shape(gtau) == [n1-a1-b1,n2-a2-b2] ) )

    ! DEBUG mode, Hermite boundary conditions:
    ! Verify that required arrays are passed and allocated with correct shape
    if (a1 > 0) then
      SLL_ASSERT(   present( boundary_data ) )
      SLL_ASSERT( allocated( boundary_data % derivs_x1_min )  )
      SLL_ASSERT(      size( boundary_data % derivs_x1_min, 1 ) == a1       )
      SLL_ASSERT(      size( boundary_data % derivs_x1_min, 2 ) == n2-a2-b2 )
    end if

    if (b1 > 0) then
      SLL_ASSERT(   present( boundary_data ) )
      SLL_ASSERT( allocated( boundary_data % derivs_x1_max )  )
      SLL_ASSERT(      size( boundary_data % derivs_x1_max, 1 ) == b1       )
      SLL_ASSERT(      size( boundary_data % derivs_x1_max, 2 ) == n2-a2-b2 )
    end if

    if (a2 > 0) then
      SLL_ASSERT(   present( boundary_data ) )
      SLL_ASSERT( allocated( boundary_data % derivs_x2_min )  )
      SLL_ASSERT(      size( boundary_data % derivs_x2_min, 1 ) == a2       )
      SLL_ASSERT(      size( boundary_data % derivs_x2_min, 2 ) == n1-a1-b1 )
    end if

    if (b2 > 0) then
      SLL_ASSERT(   present( boundary_data ) )
      SLL_ASSERT( allocated( boundary_data % derivs_x2_max )  )
      SLL_ASSERT(      size( boundary_data % derivs_x2_max, 1 ) == b2       )
      SLL_ASSERT(      size( boundary_data % derivs_x2_max, 2 ) == n1-a1-b1 )
    end if

    if (a1 > 0 .and. a2 > 0) then
      SLL_ASSERT( allocated( boundary_data % mixed_derivs_a )  )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_a, 1 ) == a1 )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_a, 2 ) == a2 )
    end if

    if (b1 > 0 .and. a2 > 0) then
      SLL_ASSERT( allocated( boundary_data % mixed_derivs_b )  )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_b, 1 ) == b1 )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_b, 2 ) == a2 )
    end if

    if (a1 > 0 .and. b2 > 0) then
      SLL_ASSERT( allocated( boundary_data % mixed_derivs_c )  )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_c, 1 ) == a1 )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_c, 2 ) == b2 )
    end if

    if (b1 > 0 .and. b2 > 0) then
      SLL_ASSERT( allocated( boundary_data % mixed_derivs_d )  )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_d, 1 ) == b1 )
      SLL_ASSERT(      size( boundary_data % mixed_derivs_d, 2 ) == b2 )
    end if

    ! Copy interpolation data onto w array
    w(1+a1:n1-b1,1+a2:n2-b2) = gtau(:,:)

    ! Hermite BCs: store boundary data in appropriate chunk of w array
    if (present( boundary_data )) then

      if (a1 > 0) w(      1:a1,1+a2:n2-b2) = boundary_data % derivs_x1_min(1:a1,:)
      if (b1 > 0) w(n1-b1+1:n1,1+a2:n2-b2) = boundary_data % derivs_x1_max(1:b1,:)

      if (a2 > 0) w(1+a1:n1-b1,      1:a2) = transpose( boundary_data % derivs_x2_min(1:a2,:) )
      if (b2 > 0) w(1+a1:n1-b1,n2-b2+1:n2) = transpose( boundary_data % derivs_x2_max(1:b2,:) )

      if (a1 > 0 .and. a2 > 0) w(      1:a1,      1:a2) = boundary_data % mixed_derivs_a(:,:)
      if (b1 > 0 .and. a2 > 0) w(n1-b1+1:n1,      1:a2) = boundary_data % mixed_derivs_b(:,:)
      if (a1 > 0 .and. b2 > 0) w(      1:a1,n2-b2+1:n2) = boundary_data % mixed_derivs_c(:,:)
      if (b1 > 0 .and. b2 > 0) w(n1-b1+1:n1,n2-b2+1:n2) = boundary_data % mixed_derivs_d(:,:)

    end if

    ! Cycle over x2 position (or order of x2-derivative at boundary)
    ! and interpolate f along x1 direction. Store coefficients in bwork
    do i2 = 1, n2

      call self%bs1%compute_interpolant( &
        gtau        = w(   1+a1:n1-b1,i2) , &
        derivs_xmin = w(      1:a1   ,i2) , &
        derivs_xmax = w(n1-b1+1:n1   ,i2) )

!      w(1:n1,i2) = self%bs1%bcoef(1:n1)
      associate (bcoef1 => self % bs1 % get_coeff())
      w(1:n1,i2) = bcoef1(1:n1)
      end associate

    end do

    wt = transpose( w )

    ! Cycle over x1 position (or order of x1-derivative at boundary)
    ! and interpolate w along x2 direction. Store coefficients in bcoef
    do i1 = 1, n1

      call self%bs2%compute_interpolant( &
        gtau        = wt(   1+a2:n2-b2,i1) , &
        derivs_xmin = wt(      1:a2   ,i1) , &
        derivs_xmax = wt(n2-b2+1:n2   ,i1) )

!      wt(1:n2,i1) = self%bs2%bcoef(1:n2)
      associate (bcoef2 => self % bs2 % get_coeff())
      wt(1:n2,i1) = bcoef2(1:n2)
      end associate

    end do

    ! x1-periodic only: "wrap around" coefficients onto extended array
    if (self%bc_xmin(1) == sll_p_periodic) then
      associate( g1 => 1+self%deg(1)/2 )
        wt(:,1-g1:0    ) = wt(:,n1-g1+1:n1)
        wt(:,n1+1:n1+g1) = wt(:,      1:g1)
      end associate
    end if

    w = transpose( wt )

    ! x2-periodic only: "wrap around" coefficients onto extended array
    if (self%bc_xmin(2) == sll_p_periodic) then
      associate( g2 => 1+self%deg(2)/2 )
        w(:,1-g2:0    ) = w(:,n2-g2+1:n2)
        w(:,n2+1:n2+g2) = w(:,      1:g2)
      end associate
    end if

    end associate

  end subroutine s_spline_2d_uniform__compute_interpolant

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_2d_uniform__eval( self, x1, x2 ) result( y )

    class(sll_t_spline_2d_uniform), intent(in) :: self
    real(wp)                      , intent(in) :: x1
    real(wp)                      , intent(in) :: x2
    real(wp) :: y

    integer  :: a1, a2
    integer  :: b1, b2
    integer  :: k1, k2
    integer  :: ic1, ic2
    real(wp) :: x1_offset, x2_offset

    ! Automatic arrays
    real(wp) :: values1(1+self%bs1%deg)
    real(wp) :: values2(1+self%bs2%deg)

    ! Find 2D cell (ic1,ic2) and normalized offset therein
    call self % bs1 % get_icell_and_offset( x1, ic1, x1_offset )
    call self % bs2 % get_icell_and_offset( x2, ic2, x2_offset )

    ! Compute arrays v1 and v2 of B-spline values
    call sll_s_uniform_bsplines_eval_basis( self%bs1%deg, x1_offset, values1 )
    call sll_s_uniform_bsplines_eval_basis( self%bs2%deg, x2_offset, values2 )

    ! Calculate (matrix) block C of B-spline coefficients to be used
    a1 = ic1 - self%bs1%offset;  b1 = ic1 - self%bs1%offset + self%bs1%deg
    a2 = ic2 - self%bs2%offset;  b2 = ic2 - self%bs2%offset + self%bs2%deg

    ! Compute scalar product <v1,v2> = (v1^T)*C*v2
    associate( bcoef => self%bcoef(a1:b1,a2:b2) )
      y = 0.0_f64
      do k2 = 1, 1+self%bs2%deg
        do k1 = 1, 1+self%bs1%deg
          y  = y + bcoef(k1,k2) * values1(k1) * values2(k2)
        end do
      end do
    end associate

  end function f_spline_2d_uniform__eval


  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_2d_uniform__eval_deriv_x1( self, x1, x2 ) result( y )

    class(sll_t_spline_2d_uniform), intent(in) :: self
    real(wp)                      , intent(in) :: x1
    real(wp)                      , intent(in) :: x2
    real(wp) :: y

    integer  :: a1, a2
    integer  :: b1, b2
    integer  :: k1, k2
    integer  :: ic1, ic2
    real(wp) :: x1_offset, x2_offset

    ! Automatic arrays
    real(wp) :: derivs1(1+self%bs1%deg)
    real(wp) :: values2(1+self%bs2%deg)

    ! Find 2D cell (ic1,ic2) and normalized offset therein
    call self % bs1 % get_icell_and_offset( x1, ic1, x1_offset )
    call self % bs2 % get_icell_and_offset( x2, ic2, x2_offset )

    ! Compute arrays d1 and v2 of B-spline values
    call sll_s_uniform_bsplines_eval_deriv( self%bs1%deg, x1_offset, derivs1 )
    call sll_s_uniform_bsplines_eval_basis( self%bs2%deg, x2_offset, values2 )

    ! Calculate (matrix) block C of B-spline coefficients to be used
    a1 = ic1 - self%bs1%offset;  b1 = ic1 - self%bs1%offset + self%bs1%deg
    a2 = ic2 - self%bs2%offset;  b2 = ic2 - self%bs2%offset + self%bs2%deg

    ! Compute scalar product <d1,v2> = (d1^T)*C*v2
    associate( bcoef => self%bcoef(a1:b1,a2:b2) )
      y = 0.0_f64
      do k2 = 1, 1+self%bs2%deg
        do k1 = 1, 1+self%bs1%deg
          y  = y + bcoef(k1,k2) * derivs1(k1) * values2(k2)
        end do
      end do
    end associate

    ! When using uniform B-splines, the cell size is normalized between 0 and 1,
    ! hence the derivative must be divided by dx(1) to scale back the result.
    y = y * self%inv_dx(1)

  end function f_spline_2d_uniform__eval_deriv_x1

  !-----------------------------------------------------------------------------
  SLL_PURE function f_spline_2d_uniform__eval_deriv_x2( self, x1, x2 ) result( y )

    class(sll_t_spline_2d_uniform), intent(in) :: self
    real(wp)                      , intent(in) :: x1
    real(wp)                      , intent(in) :: x2
    real(wp) :: y

    integer  :: a1, a2
    integer  :: b1, b2
    integer  :: k1, k2
    integer  :: ic1, ic2
    real(wp) :: x1_offset, x2_offset

    ! Automatic arrays
    real(wp) :: values1(1+self%bs1%deg)
    real(wp) :: derivs2(1+self%bs2%deg)

    ! Find 2D cell (ic1,ic2) and normalized offset therein
    call self % bs1 % get_icell_and_offset( x1, ic1, x1_offset )
    call self % bs2 % get_icell_and_offset( x2, ic2, x2_offset )

    ! Compute arrays v1 and d2 of B-spline values
    call sll_s_uniform_bsplines_eval_basis( self%bs1%deg, x1_offset, values1 )
    call sll_s_uniform_bsplines_eval_deriv( self%bs2%deg, x2_offset, derivs2 )

    ! Calculate (matrix) block C of B-spline coefficients to be used
    a1 = ic1 - self%bs1%offset;  b1 = ic1 - self%bs1%offset + self%bs1%deg
    a2 = ic2 - self%bs2%offset;  b2 = ic2 - self%bs2%offset + self%bs2%deg

    ! Compute scalar product <v1,d2> = (v1^T)*C*d2
    associate( bcoef => self%bcoef(a1:b1,a2:b2) )
      y = 0.0_f64
      do k2 = 1, 1+self%bs2%deg
        do k1 = 1, 1+self%bs1%deg
          y  = y + bcoef(k1,k2) * values1(k1) * derivs2(k2)
        end do
      end do
    end associate

    ! When using uniform B-splines, the cell size is normalized between 0 and 1,
    ! hence the derivative must be divided by dx(2) to scale back the result.
    y = y * self%inv_dx(2)

  end function f_spline_2d_uniform__eval_deriv_x2

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_2d_uniform__eval_array( self, x1, x2, y )

    class(sll_t_spline_2d_uniform), intent(in   ) :: self
    real(wp)                      , intent(in   ) :: x1(:,:)
    real(wp)                      , intent(in   ) :: x2(:,:)
    real(wp)                      , intent(  out) :: y (:,:)

    integer :: i1, i2
    integer :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = f_spline_2d_uniform__eval( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine s_spline_2d_uniform__eval_array

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_2d_uniform__eval_array_deriv_x1( self, x1, x2, y )

    class(sll_t_spline_2d_uniform), intent(in   ) :: self
    real(wp)                      , intent(in   ) :: x1(:,:)
    real(wp)                      , intent(in   ) :: x2(:,:)
    real(wp)                      , intent(  out) :: y (:,:)

    integer :: i1, i2
    integer :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = f_spline_2d_uniform__eval_deriv_x1( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine s_spline_2d_uniform__eval_array_deriv_x1

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_spline_2d_uniform__eval_array_deriv_x2( self, x1, x2, y )

    class(sll_t_spline_2d_uniform), intent(in   ) :: self
    real(wp)                      , intent(in   ) :: x1(:,:)
    real(wp)                      , intent(in   ) :: x2(:,:)
    real(wp)                      , intent(  out) :: y (:,:)

    integer :: i1, i2
    integer :: n1, n2

    SLL_ASSERT( all( shape(x1) == shape(x2) ) )

    n1 = size(x1,1)
    n2 = size(x1,2)

    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = f_spline_2d_uniform__eval_deriv_x2( self, x1(i1,i2), x2(i1,i2) )
      end do
    end do

  end subroutine s_spline_2d_uniform__eval_array_deriv_x2

  !-----------------------------------------------------------------------------
  subroutine s_spline_2d_uniform__free( self )

    class(sll_t_spline_2d_uniform), intent(inout) :: self

    ! Deallocate local 2D arrays
    deallocate( self%bwork )
    deallocate( self%bcoef )

    ! Free memory of 1D B-splines
    call self%bs1%free()
    call self%bs2%free()

  end subroutine s_spline_2d_uniform__free

  !-----------------------------------------------------------------------------
  subroutine s_spline_2d_uniform__get_interp_points( self, tau1, tau2 )

    class(sll_t_spline_2d_uniform), intent(in   ) :: self
    real(wp),          allocatable, intent(  out) :: tau1(:)
    real(wp),          allocatable, intent(  out) :: tau2(:)

    call self % bs1 % get_interp_points( tau1 )
    call self % bs2 % get_interp_points( tau2 )

  end subroutine s_spline_2d_uniform__get_interp_points

end module sll_m_spline_2d_uniform
