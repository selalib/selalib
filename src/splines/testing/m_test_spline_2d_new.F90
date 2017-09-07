module m_test_spline_2d_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: &
    f64

  use m_analytical_profiles_2d_base, only: &
    t_profile_2d_info, &
    c_analytical_profile_2d

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite, &
    sll_p_greville

  use sll_m_bsplines_base, only: &
    sll_c_bsplines

  use sll_m_bsplines_uniform, only: &
    sll_t_bsplines_uniform

  use sll_m_bsplines_non_uniform, only: &
    sll_t_bsplines_non_uniform

  use sll_m_spline_2d_new, only: &
    sll_t_spline_2d

  use sll_m_spline_interpolator_2d, only: &
    sll_t_spline_interpolator_2d, &
    sll_t_spline_2d_boundary_data

  use sll_m_timer, only: &
    sll_t_time_mark, &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_between

  implicit none

  public :: &
    t_spline_2d_test_facility

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Type for running test
  type :: t_spline_2d_test_facility

    class(c_analytical_profile_2d), pointer :: profile_2d

    class(sll_c_bsplines), allocatable :: bsplines_x1
    class(sll_c_bsplines), allocatable :: bsplines_x2
    type(sll_t_spline_2d)              :: spline_2d
    type(sll_t_spline_interpolator_2d) :: spline_interpolator_2d

    real(wp), allocatable ::  tau1(:)   ! interp. points, x1 coord
    real(wp), allocatable ::  tau2(:)   ! interp. points, x2 coord
    real(wp), allocatable :: gtau (:,:) ! profile values at tau

    real(wp) :: time_init
    real(wp) :: time_compute_interpolant
    real(wp) :: time_eval
    real(wp) :: time_eval_array
    real(wp) :: time_eval_diff1
    real(wp) :: time_eval_diff1_array
    real(wp) :: time_eval_diff2
    real(wp) :: time_eval_diff2_array

  contains

    procedure :: init
    procedure :: free
    procedure :: check_equivalence_scalar_array_methods
    procedure :: evaluate_on_2d_grid
    procedure :: evaluate_at_interpolation_points
    procedure :: evaluate_grad_on_2d_grid

  end type t_spline_2d_test_facility

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  subroutine init( self, profile_2d, degree, ncells, bc_xmin, bc_xmax, breaks1, breaks2 )

    class(t_spline_2d_test_facility), intent(  out)         :: self
    class(c_analytical_profile_2d  ), intent(in   ), target :: profile_2d
    integer                         , intent(in   )         :: degree (2)
    integer                         , intent(in   )         :: ncells (2)
    integer                         , intent(in   )         :: bc_xmin(2)
    integer                         , intent(in   )         :: bc_xmax(2)
    real(wp),               optional, intent(in   )         :: breaks1(:)
    real(wp),               optional, intent(in   )         :: breaks2(:)

    type(t_profile_2d_info) :: info
    integer                 :: i1, j1
    integer                 :: i2, j2
    real(wp)                :: dx1, dx2
    logical                 :: periodic(2)

    type(sll_t_spline_2d_boundary_data) :: boundary_data

    type(sll_t_time_mark) :: t0, t1, t2, t3

    ! Store pointer to 2D profile
    self % profile_2d => profile_2d

    ! Extract information about 2D analytical profile
    call self % profile_2d % get_info( info )

    ! TODO: breaks1 and breaks2 should be read from input arguments
    dx1 = (info%x1_max-info%x1_min)/real(ncells(1),wp)
    dx2 = (info%x2_max-info%x2_min)/real(ncells(2),wp)

    call sll_s_set_time_mark( t0 )

    !TODO: add proper checks
    if (bc_xmin(1) == sll_p_periodic) then
      periodic(1) = .true.
    else
      periodic(1) = .false.
    end if

    !TODO: add proper checks
    if (bc_xmin(2) == sll_p_periodic) then
      periodic(2) = .true.
    else
      periodic(2) = .false.
    end if

    ! Allocate B-splines (uniform or non-uniform depending on input)
    if (present( breaks1 ) .and. size( breaks1 )>0 ) then
      allocate( sll_t_bsplines_non_uniform :: self % bsplines_x1 )
    else 
      allocate( sll_t_bsplines_uniform     :: self % bsplines_x1 )
    end if

    ! Allocate B-splines (uniform or non-uniform depending on input)
    if (present( breaks2 ) .and. size( breaks2 )>0 ) then
      allocate( sll_t_bsplines_non_uniform :: self % bsplines_x2 )
    else 
      allocate( sll_t_bsplines_uniform     :: self % bsplines_x2 )
    end if

    ! Initialize B-splines (uniform or non-uniform depending on allocated type)
    select type( bsplines_x1 => self % bsplines_x1 )
      type is( sll_t_bsplines_non_uniform )
        call bsplines_x1 % init( degree(1), periodic(1), breaks1 )
      type is( sll_t_bsplines_uniform     )
        call bsplines_x1 % init( degree(1), periodic(1), info%x1_min, info%x1_max, ncells(1) )
    end select

    ! Initialize B-splines (uniform or non-uniform depending on allocated type)
    select type( bsplines_x2 => self % bsplines_x2 )
      type is( sll_t_bsplines_non_uniform )
        call bsplines_x2 % init( degree(2), periodic(2), breaks2 )
      type is( sll_t_bsplines_uniform     )
        call bsplines_x2 % init( degree(2), periodic(2), info%x2_min, info%x2_max, ncells(2) )
    end select

    ! Initialize 2D spline
    call self % spline_2d % init( self % bsplines_x1, self % bsplines_x2 )

    ! Initialize 2D interpolator
    call self % spline_interpolator_2d % init( &
      self % bsplines_x1, &
      self % bsplines_x2, &
      bc_xmin(:)        , &
      bc_xmax(:) )

    call sll_s_set_time_mark( t1 )

    ! Get x1 and x2 coordinates of interpolation points
    call self % spline_interpolator_2d % get_interp_points( self % tau1, self % tau2 )

    associate( tau1 => self%tau1, tau2 => self%tau2 )

    associate( nipts1 => size( tau1 )         , & ! number of interpolation pts
               nipts2 => size( tau2 )         , & ! number of interpolation pts
               s1     => 1-modulo(degree(1),2), & ! shift = 1 for even order, 0 for odd order
               s2     => 1-modulo(degree(2),2), & ! shift = 1 for even order, 0 for odd order
               dh1    => degree(1)/2          , &
               dh2    => degree(2)/2          )

    ! Evaluate analytical profile at interpolation points
    allocate( self % gtau (nipts1,nipts2) )
    do i2 = 1, nipts2
      do i1 = 1, nipts1
        self % gtau(i1,i2) = self % profile_2d % eval( tau1(i1), tau2(i2) )
      end do
    end do

    ! If needed, evaluate x1 derivatives at (x1_min,x2)
    if (bc_xmin(1) == sll_p_hermite) then
      allocate( boundary_data % derivs_x1_min (dh1, nipts2) )
      do i2 = 1, nipts2
        do j1 = 1, dh1
          boundary_data % derivs_x1_min(j1,i2) = &
            self % profile_2d % eval( info%x1_min, tau2(i2), diff_x1=j1-s1 )
        end do
      end do
    end if

    ! If needed, evaluate x1 derivatives at (x1_max,x2)
    if (bc_xmax(1) == sll_p_hermite) then
      allocate( boundary_data % derivs_x1_max (dh1, nipts2) )
      do i2 = 1, nipts2
        do j1 = 1, dh1
          boundary_data % derivs_x1_max(j1,i2) = &
            self % profile_2d % eval( info%x1_max, tau2(i2), diff_x1=j1-s1 )
        end do
      end do
    end if

    ! If needed, evaluate x2 derivatives at (x1,x2_min)
    if (bc_xmin(2) == sll_p_hermite) then
      allocate( boundary_data % derivs_x2_min (dh2, nipts1) )
      do i1 = 1, nipts1
        do j2 = 1, dh2
          boundary_data % derivs_x2_min(j2,i1) = &
            self % profile_2d % eval( tau1(i1), info%x2_min, diff_x2=j2-s2 )
        end do
      end do
    end if

    ! If needed, evaluate x2 derivatives at (x1,x2_max)
    if (bc_xmax(2) == sll_p_hermite) then
      allocate( boundary_data % derivs_x2_max (dh2, nipts1) )
      do i1 = 1, nipts1
        do j2 = 1, dh2
          boundary_data % derivs_x2_max(j2,i1) = &
            self % profile_2d % eval( tau1(i1), info%x2_max, diff_x2=j2-s2 )
        end do
      end do
    end if

    ! If needed, evaluate (x1,x2) mixed derivatives at corner A
    if (all( [bc_xmin(1),bc_xmin(2)] == sll_p_hermite )) then
      allocate( boundary_data % mixed_derivs_a (dh1,dh2) )
      do j2 = 1, dh2
        do j1 = 1, dh1
          boundary_data % mixed_derivs_a(j1,j2) = &
            profile_2d % eval( info%x1_min, info%x2_min, diff_x1=j1-s1, diff_x2=j2-s2 )
        end do
      end do
    end if

    ! If needed, evaluate (x1,x2) mixed derivatives at corner B
    if (all( [bc_xmax(1),bc_xmin(2)] == sll_p_hermite )) then
      allocate( boundary_data % mixed_derivs_b (dh1,dh2) )
      do j2 = 1, dh2
        do j1 = 1, dh1
          boundary_data % mixed_derivs_b(j1,j2) = &
            profile_2d % eval( info%x1_max, info%x2_min, diff_x1=j1-s1, diff_x2=j2-s2 )
        end do
      end do
    end if

    ! If needed, evaluate (x1,x2) mixed derivatives at corner C
    if (all( [bc_xmin(1),bc_xmax(2)] == sll_p_hermite )) then
      allocate( boundary_data % mixed_derivs_c (dh1,dh2) )
      do j2 = 1, dh2
        do j1 = 1, dh1
          boundary_data % mixed_derivs_c(j1,j2) = &
            profile_2d % eval( info%x1_min, info%x2_max, diff_x1=j1-s1, diff_x2=j2-s2 )
        end do
      end do
    end if

    ! If needed, evaluate (x1,x2) mixed derivatives at corner D
    if (all( [bc_xmax(1),bc_xmax(2)] == sll_p_hermite )) then
      allocate( boundary_data % mixed_derivs_d (dh1,dh2) )
      do j2 = 1, dh2
        do j1 = 1, dh1
          boundary_data % mixed_derivs_d(j1,j2) = &
            profile_2d % eval( info%x1_max, info%x2_max, diff_x1=j1-s1, diff_x2=j2-s2 )
        end do
      end do
    end if

    end associate ! nipts1, nipts2, s1, s2, dh1, dh2
    end associate ! tau1, tau2

    ! Compute 2D spline that interpolates analytical 2D profile at points above
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call sll_s_set_time_mark( t2 )

    call self % spline_interpolator_2d % compute_interpolant( &
      self % spline_2d, &
      self % gtau     , &
      boundary_data )

    call sll_s_set_time_mark( t3 )

    ! Timings (set to -1 values not yet available)
    self % time_init                =  sll_f_time_elapsed_between( t0, t1 )
    self % time_compute_interpolant =  sll_f_time_elapsed_between( t2, t3 )
    self % time_eval                = -1.0_wp
    self % time_eval_array          = -1.0_wp
    self % time_eval_diff1          = -1.0_wp
    self % time_eval_diff1_array    = -1.0_wp
    self % time_eval_diff2          = -1.0_wp
    self % time_eval_diff2_array    = -1.0_wp

  end subroutine init

  !-----------------------------------------------------------------------------
  subroutine check_equivalence_scalar_array_methods( self, equiv )

    class(t_spline_2d_test_facility), intent(inout) :: self
    logical                         , intent(  out) :: equiv(3)

    type(t_profile_2d_info) :: info
    real(wp), allocatable   :: x1(:,:)
    real(wp), allocatable   :: x2(:,:)
    real(wp), allocatable   :: y (:,:)
    real(wp), allocatable   :: ya(:,:)
    real(wp)                :: r
    integer                 :: i1, i2
    type(sll_t_time_mark)   :: t0, t1, t2

    ! Choose size of 2D grid
    integer , parameter :: n1   = 77
    integer , parameter :: n2   = 90
    real(wp), parameter :: npts = real( n1*n2, wp )

    ! Allocate local arrays
    allocate( x1(n1,n2) )
    allocate( x2(n1,n2) )
    allocate( y (n1,n2) )
    allocate( ya(n1,n2) )

    ! Get info from profile
    call self % profile_2d % get_info( info )

    ! Generate random grids in x1 and x2
    do i2 = 1, n2
      do i1 = 1, n1
        call random_number( r )
        x1(i1,i2) = info%x1_min*(1.0_wp-r) + info%x1_max*r
        x2(i1,i2) = info%x2_min*(1.0_wp-r) + info%x2_max*r
      end do
    end do

    ! Compare:
    !   . sll_f_spline_2d_non_uniform_eval
    !   . sll_s_spline_2d_non_uniform_eval_array
    call sll_s_set_time_mark( t0 )
    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = self % spline_2d % eval( x1(i1,i2), x2(i1,i2) )
      end do
    end do
    call sll_s_set_time_mark( t1 )
    call self % spline_2d % eval_array( x1, x2, ya )
    call sll_s_set_time_mark( t2 )

    self%time_eval       = sll_f_time_elapsed_between( t0, t1 ) / npts
    self%time_eval_array = sll_f_time_elapsed_between( t1, t2 ) / npts
    equiv(1) = all( y == ya )

    ! Compare:
    !   . sll_f_spline_2d_non_uniform_eval_deriv_x1
    !   . sll_s_spline_2d_non_uniform_eval_array_deriv_x1
    call sll_s_set_time_mark( t0 )
    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = self % spline_2d % eval_deriv_x1( x1(i1,i2), x2(i1,i2) )
      end do
    end do
    call sll_s_set_time_mark( t1 )
    call self % spline_2d % eval_array_deriv_x1( x1, x2, ya )
    call sll_s_set_time_mark( t2 )

    self%time_eval_diff1       = sll_f_time_elapsed_between( t0, t1 ) / npts
    self%time_eval_diff1_array = sll_f_time_elapsed_between( t1, t2 ) / npts
    equiv(2) = all( y == ya )

    ! Compare:
    !   . sll_f_spline_2d_non_uniform_eval_deriv_x2
    !   . sll_s_spline_2d_non_uniform_eval_array_deriv_x2
    call sll_s_set_time_mark( t0 )
    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = self % spline_2d % eval_deriv_x2( x1(i1,i2), x2(i1,i2) )
      end do
    end do
    call sll_s_set_time_mark( t1 )
    call self % spline_2d % eval_array_deriv_x2( x1, x2, ya )
    call sll_s_set_time_mark( t2 )

    self%time_eval_diff2       = sll_f_time_elapsed_between( t0, t1 ) / npts
    self%time_eval_diff2_array = sll_f_time_elapsed_between( t1, t2 ) / npts
    equiv(3) = all( y == ya )

    ! Deallocate local arrays
    deallocate( x1 )
    deallocate( x2 )
    deallocate( y  )
    deallocate( ya )

  end subroutine check_equivalence_scalar_array_methods

  !-----------------------------------------------------------------------------
  subroutine evaluate_at_interpolation_points( self, max_norm_error )

    class(t_spline_2d_test_facility), intent(in   ) :: self
    real(wp)                        , intent(  out) :: max_norm_error

    integer  :: i1, i2
    real(wp) :: error

    ! Evaluate 2D spline at interpolation points:
    ! interpolation values should be obtained
    max_norm_error = 0.0_wp
    do i2 = 1, size( self%tau2 )
      do i1 = 1, size( self%tau1 )
        error = self % gtau(i1,i2) &
              - self % spline_2d % eval( self%tau1(i1), self%tau2(i2) )
        max_norm_error = max( max_norm_error, abs( error ) )
      end do
    end do

  end subroutine evaluate_at_interpolation_points

  !-----------------------------------------------------------------------------
  subroutine evaluate_on_2d_grid( self, x1, x2, max_norm_error )

    class(t_spline_2d_test_facility), intent(inout) :: self
    real(wp)                        , intent(in   ) :: x1(:,:)
    real(wp)                        , intent(in   ) :: x2(:,:)
    real(wp)                        , intent(  out) :: max_norm_error

    integer               :: i1, i2
    integer               :: n1, n2
    real(wp), allocatable :: y(:,:)
    real(wp)              :: error
    type(sll_t_time_mark) :: t0, t1

    SLL_ASSERT( all( shape( x1 ) == shape( x2 ) ) )

    n1 = size( x1, 1 )
    n2 = size( x2, 2 )

    ! Array of spline values
    allocate( y(n1,n2) )

    call sll_s_set_time_mark( t0 )

    ! Evaluate 2D spline on given grid
    call self % spline_2d % eval_array( x1, x2, y )

    call sll_s_set_time_mark( t1 )
    self % time_eval_array = sll_f_time_elapsed_between(t0,t1)/real(n1*n2,wp)

    ! Compare spline values to analytical profile and compute max norm of error
    max_norm_error = 0.0_wp
    do i2 = 1, n2
      do i1 = 1, n1
        error = self % profile_2d % eval( x1(i1,i2), x2(i1,i2) ) - y(i1,i2)
        max_norm_error = max( max_norm_error, abs( error ) )
      end do
    end do

    ! Free local memory
    deallocate( y )

  end subroutine evaluate_on_2d_grid

  !-----------------------------------------------------------------------------
  subroutine evaluate_grad_on_2d_grid( self, &
      x1, &
      x2, &
      max_norm_error_diff_x1, &
      max_norm_error_diff_x2 )

    class(t_spline_2d_test_facility), intent(inout) :: self
    real(wp)                        , intent(in   ) :: x1(:,:)
    real(wp)                        , intent(in   ) :: x2(:,:)
    real(wp)                        , intent(  out) :: max_norm_error_diff_x1
    real(wp)                        , intent(  out) :: max_norm_error_diff_x2

    integer               :: i1, i2
    integer               :: n1, n2
    real(wp), allocatable :: y(:,:)
    real(wp)              :: error
    type(sll_t_time_mark) :: t0, t1

    SLL_ASSERT( all( shape( x1 ) == shape( x2 ) ) )

    n1 = size( x1, 1 )
    n2 = size( x2, 2 )

    ! Array of spline values
    allocate( y(n1,n2) )

    ! x1-derivative: compare spline to analytical profile
    !----------------------------------------------------
    call sll_s_set_time_mark( t0 )

    call self % spline_2d % eval_array_deriv_x1( x1, x2, y )

    call sll_s_set_time_mark( t1 )
    self % time_eval_diff1_array = sll_f_time_elapsed_between(t0,t1)/real(n1*n2,wp)

    max_norm_error_diff_x1 = 0.0_wp
    do i2 = 1, n2
      do i1 = 1, n1
        error = self % profile_2d % eval( x1(i1,i2), x2(i1,i2), diff_x1=1 ) - y(i1,i2)
        max_norm_error_diff_x1 = max( max_norm_error_diff_x1, abs( error ) )
      end do
    end do

    ! x2-derivative: compare spline to analytical profile
    !----------------------------------------------------
    call sll_s_set_time_mark( t0 )

    call self % spline_2d % eval_array_deriv_x2( x1, x2, y )

    call sll_s_set_time_mark( t1 )
    self % time_eval_diff2_array = sll_f_time_elapsed_between(t0,t1)/real(n1*n2,wp)

    max_norm_error_diff_x2 = 0.0_wp
    do i2 = 1, n2
      do i1 = 1, n1
        error = self % profile_2d % eval( x1(i1,i2), x2(i1,i2), diff_x2=1 ) - y(i1,i2)
        max_norm_error_diff_x2 = max( max_norm_error_diff_x2, abs( error ) )
      end do
    end do

    ! Free local memory
    deallocate( y )

  end subroutine evaluate_grad_on_2d_grid

  !-----------------------------------------------------------------------------
  subroutine free( self )

    class(t_spline_2d_test_facility), intent(inout) :: self

    ! Free spline memory
    call self % spline_2d % free()

  end subroutine free

end module m_test_spline_2d_new
