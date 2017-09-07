module sll_m_spline_interpolator_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision     , only: f64
  use sll_m_bsplines_base         , only: sll_c_bsplines
  use sll_m_spline_1d_new         , only: sll_t_spline_1d
  use sll_m_spline_2d_new         , only: sll_t_spline_2d
  use sll_m_spline_interpolator_1d, only: sll_t_spline_interpolator_1d

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite, &
    sll_p_greville

  implicit none

  public :: &
    sll_t_spline_interpolator_2d, &
    sll_t_spline_2d_boundary_data

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64


  !> Container for boundary condition data
  !>
  !>  x2_max  ____________
  !>         |            |
  !>         | c        d |
  !>         |            |
  !>         |            |
  !>         |            |
  !>         | a        b |
  !>  x2_min |____________|
  !>       x1_min       x1_max
  !>
  type :: sll_t_spline_2d_boundary_data
    real(wp), allocatable :: derivs_x1_min (:,:)
    real(wp), allocatable :: derivs_x1_max (:,:)
    real(wp), allocatable :: derivs_x2_min (:,:)
    real(wp), allocatable :: derivs_x2_max (:,:)
    real(wp), allocatable :: mixed_derivs_a(:,:)
    real(wp), allocatable :: mixed_derivs_b(:,:)
    real(wp), allocatable :: mixed_derivs_c(:,:)
    real(wp), allocatable :: mixed_derivs_d(:,:)
  end type sll_t_spline_2d_boundary_data


  !> 2D spline interpolator
  type :: sll_t_spline_interpolator_2d

    ! Private attributes
    class(sll_c_bsplines),     pointer, private :: bspl1 => null()
    class(sll_c_bsplines),     pointer, private :: bspl2 => null()
    integer                           , private ::  bc_xmin(2)
    integer                           , private ::  bc_xmax(2)
    integer                           , private :: nbc_xmin(2)
    integer                           , private :: nbc_xmax(2)
    type(sll_t_spline_1d)             , private :: spline1
    type(sll_t_spline_1d)             , private :: spline2
    type(sll_t_spline_interpolator_1d), private :: interp1
    type(sll_t_spline_interpolator_1d), private :: interp2
    real(wp),              allocatable, private :: bwork(:,:)

  contains

    procedure :: init                => s_spline_interpolator_2d__init
    procedure :: free                => s_spline_interpolator_2d__free
    procedure :: get_interp_points   => s_spline_interpolator_2d__get_interp_points
    procedure :: compute_interpolant => s_spline_interpolator_2d__compute_interpolant

  end type sll_t_spline_interpolator_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief Initialize a 2D spline interpolation object
  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_2d__init( &
    self   , &
    bspl1  , &
    bspl2  , &
    bc_xmin, &
    bc_xmax )

    class(sll_t_spline_interpolator_2d), intent(  out) :: self
    class(sll_c_bsplines),       target, intent(in   ) :: bspl1
    class(sll_c_bsplines),       target, intent(in   ) :: bspl2
    integer                            , intent(in   ) :: bc_xmin(2)
    integer                            , intent(in   ) :: bc_xmax(2)

    integer :: g1, g2

    ! Save pointers to B-splines
    ! (later needed to verify 2D spline input to 'compute_interpolant')
    self % bspl1 => bspl1
    self % bspl2 => bspl2

    ! Initialize 1D splines along x1 and x2
    call self % spline1 % init( bspl1 )
    call self % spline2 % init( bspl2 )

    ! Initialize 1D interpolators along x1 and x2
    call self % interp1 % init( bspl1, bc_xmin(1), bc_xmax(1) )
    call self % interp2 % init( bspl2, bc_xmin(2), bc_xmax(2) )

    associate( n1     =>  bspl1 % nbasis, &
               n2     =>  bspl2 % nbasis, &
               degree => [bspl1 % degree, bspl2 % degree] )

      g1 = merge( 1+degree(1)/2, 0, bc_xmin(1) == sll_p_periodic )
      g2 = merge( 1+degree(2)/2, 0, bc_xmin(2) == sll_p_periodic )

!      ! Save data into type
!      self%deg     = degree
!      self%n       = [n1,n2]
      self%bc_xmin = bc_xmin
      self%bc_xmax = bc_xmax

      ! Allocate work array for interpolation.
      ! In case of periodic BCs, a larger array of coefficients is used in order
      ! to avoid a loop with calls to the "mod( , )" function at evaluation.
      allocate( self%bwork(1-g2:n2+g2,1-g1:n1+g1) )
      self%bwork = 0.0_f64

      ! Calculate number of additional boundary data on each side of domain
      ! (i.e. derivatives for Hermite BC)
      self%nbc_xmin = merge( degree/2, 0, bc_xmin == sll_p_hermite )
      self%nbc_xmax = merge( degree/2, 0, bc_xmax == sll_p_hermite )

    end associate

  end subroutine s_spline_interpolator_2d__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_2d__free( self )

    class(sll_t_spline_interpolator_2d), intent(inout) :: self

    ! Deallocate local 2D arrays
    deallocate( self % bwork )

    ! Free memory of 1D B-splines
    call self % bspl1 % free()
    call self % bspl2 % free()

  end subroutine s_spline_interpolator_2d__free

  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_2d__get_interp_points( self, tau1, tau2 )

    class(sll_t_spline_interpolator_2d), intent(in   ) :: self
    real(wp),               allocatable, intent(  out) :: tau1(:)
    real(wp),               allocatable, intent(  out) :: tau2(:)

    call self % interp1 % get_interp_points( tau1 )
    call self % interp2 % get_interp_points( tau2 )

  end subroutine s_spline_interpolator_2d__get_interp_points

  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_2d__compute_interpolant( self, &
      spline, gtau, boundary_data )

    class(sll_t_spline_interpolator_2d), intent(inout)           :: self
    type (sll_t_spline_2d)             , intent(inout)           :: spline
    real(wp)                           , intent(in   )           :: gtau(:,:)
    type(sll_t_spline_2d_boundary_data), intent(in   ), optional :: boundary_data

    character(len=*), parameter :: &
      this_sub_name = "sll_t_spline_interpolator_2d % compute_interpolant"

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

    associate( w  => spline%bcoef     , &
               wt => self%bwork       , &
               n1 => self%bspl1%nbasis, &
               n2 => self%bspl2%nbasis, &
               a1 => self%nbc_xmin(1) , &
               a2 => self%nbc_xmin(2) , &
               b1 => self%nbc_xmax(1) , &
               b2 => self%nbc_xmax(2) )

    SLL_ASSERT( all( shape(gtau) == [n1-a1-b1,n2-a2-b2] ) )
    SLL_ASSERT( spline % belongs_to_space( self % bspl1, self % bspl2 ) )

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

      call self % interp1 % compute_interpolant( &
        spline      = self % spline1      , &
        gtau        = w(   1+a1:n1-b1,i2) , &
        derivs_xmin = w(      1:a1   ,i2) , &
        derivs_xmax = w(n1-b1+1:n1   ,i2) )

      w(1:n1,i2) = self % spline1 % bcoef(1:n1)

    end do

    wt = transpose( w )

    ! Cycle over x1 position (or order of x1-derivative at boundary)
    ! and interpolate w along x2 direction. Store coefficients in bcoef
    do i1 = 1, n1

      call self % interp2 % compute_interpolant( &
        spline      = self % spline2       , &
        gtau        = wt(   1+a2:n2-b2,i1) , &
        derivs_xmin = wt(      1:a2   ,i1) , &
        derivs_xmax = wt(n2-b2+1:n2   ,i1) )

      wt(1:n2,i1) = self % spline2 % bcoef(1:n2)

    end do

    ! x1-periodic only: "wrap around" coefficients onto extended array
    if (self%bc_xmin(1) == sll_p_periodic) then
      associate( g1 => 1+self%bspl1%degree/2 )
        wt(:,1-g1:0    ) = wt(:,n1-g1+1:n1)
        wt(:,n1+1:n1+g1) = wt(:,      1:g1)
      end associate
    end if

    w = transpose( wt )

    ! x2-periodic only: "wrap around" coefficients onto extended array
    if (self%bc_xmin(2) == sll_p_periodic) then
      associate( g2 => 1+self%bspl2%degree/2 )
        w(:,1-g2:0    ) = w(:,n2-g2+1:n2)
        w(:,n2+1:n2+g2) = w(:,      1:g2)
      end associate
    end if

    end associate

  end subroutine s_spline_interpolator_2d__compute_interpolant


end module sll_m_spline_interpolator_2d
