module m_test_bsplines_2d
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

  use sll_m_bspline_2d, only:             &
    sll_t_bspline_2d,                     &
    sll_s_bspline_2d_init,                &
    sll_s_bspline_2d_free,                &
    sll_s_bspline_2d_compute_interpolant, &
    sll_f_bspline_2d_eval,                &
    sll_f_bspline_2d_eval_deriv_x1,       &
    sll_f_bspline_2d_eval_deriv_x2,       &
    sll_s_bspline_2d_eval_array,          &
    sll_s_bspline_2d_eval_array_deriv_x1, &
    sll_s_bspline_2d_eval_array_deriv_x2

  use sll_m_timer, only: &
    sll_t_time_mark, &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_between

  implicit none

  public :: &
    t_bspline_2d_test_facility

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Type for running test
  type :: t_bspline_2d_test_facility

    class(c_analytical_profile_2d), pointer :: profile_2d
    integer                                 :: nx1
    integer                                 :: nx2
    integer                                 :: deg1
    integer                                 :: deg2
    integer                                 :: bc1
    integer                                 :: bc2

    type(sll_t_bspline_2d)  :: bspline_2d
    real(wp), allocatable   :: gtau(:,:)  ! Profile values at interp. points

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

  end type t_bspline_2d_test_facility

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  subroutine init( self, profile_2d, nx1, nx2, deg1, deg2, bc1, bc2 )

    class(t_bspline_2d_test_facility), intent(  out)         :: self
    class(c_analytical_profile_2d   ), intent(in   ), target :: profile_2d
    integer                          , intent(in   )         :: nx1
    integer                          , intent(in   )         :: nx2
    integer                          , intent(in   )         :: deg1
    integer                          , intent(in   )         :: deg2
    integer                          , intent(in   )         :: bc1
    integer                          , intent(in   )         :: bc2

    type(t_profile_2d_info) :: info
    integer                 :: nipts1
    integer                 :: nipts2
    integer                 :: i1, j1, s1
    integer                 :: i2, j2, s2

    real(wp), allocatable :: val1_min(:,:)
    real(wp), allocatable :: val1_max(:,:)
    real(wp), allocatable :: val2_min(:,:)
    real(wp), allocatable :: val2_max(:,:)
    real(wp), allocatable :: val_corners(:,:,:)

    type(sll_t_time_mark) :: t0, t1, t2, t3

    ! Store pointer to 2D profile and input numerical parameters to spline
    self % profile_2d => profile_2d
    self % nx1        =  nx1
    self % nx2        =  nx2
    self % deg1       =  deg1
    self % deg2       =  deg2
    self % bc1        =  bc1
    self % bc2        =  bc2

    ! Extract information about 2D analytical profile
    call self % profile_2d % get_info( info )

    call sll_s_set_time_mark( t0 )

    ! Initialize 2D spline
    call sll_s_bspline_2d_init( &
      self % bspline_2d, &
      self % nx1, &
      self % nx2, &
      self % deg1, &
      self % deg2, &
      info % x1_min, &
      info % x2_min, &
      info % x1_max, &
      info % x2_max, &
      self % bc1, &  ! BC @ x1_min
      self % bc2, &  ! BC @ x2_min
      self % bc1, &  ! BC @ x1_max
      self % bc2 )   ! BC @ x2_max

    call sll_s_set_time_mark( t1 )

    ! Get spline interpolation points
    associate( tau1 => self % bspline_2d % bs1 % tau, &
               tau2 => self % bspline_2d % bs2 % tau )

    ! Store number of interpolation points
    nipts1 = size( tau1 )
    nipts2 = size( tau2 )

    ! Evaluate analytical profile at interpolation points
    allocate( self % gtau (nipts1,nipts2) )
    do i2 = 1, nipts2
      do i1 = 1, nipts1
        self % gtau(i1,i2) = self % profile_2d % eval( tau1(i1), tau2(i2) )
      end do
    end do

    ! If needed, evaluate x1 derivatives at (x1_min,x2) and (x1_max,x2)
    if (self%bc1 == sll_p_hermite) then
      allocate( val1_min (deg1/2, nipts2) )
      allocate( val1_max (deg1/2, nipts2) )
      s1 = 1-modulo(deg1,2) ! shift = 1 for even order, 0 for odd order
      do i2 = 1, nipts2
        do j1 = 1, deg1/2
          val1_min(j1,i2) = self % profile_2d % eval( info%x1_min, tau2(i2), diff_x1=j1-s1 )
          val1_max(j1,i2) = self % profile_2d % eval( info%x1_max, tau2(i2), diff_x1=j1-s1 )
        end do
      end do
    end if

    ! If needed, evaluate x2 derivatives at (x1,x2_min) and (x1,x2_max)
    if (self%bc2 == sll_p_hermite) then
      allocate( val2_min (deg2/2, nipts1) )
      allocate( val2_max (deg2/2, nipts1) )
      s2 = 1-modulo(deg2,2) ! shift = 1 for even order, 0 for odd order
      do i1 = 1, nipts1
        do j2 = 1, deg2/2
          val2_min(j2,i1) = self % profile_2d % eval( tau1(i1), info%x2_min, diff_x2=j2-s2 )
          val2_max(j2,i1) = self % profile_2d % eval( tau1(i1), info%x2_max, diff_x2=j2-s2 )
        end do
      end do
    end if

    ! If needed, evaluate (x1,x2) mixed derivatives at 4 corners
    if (self%bc1 == sll_p_hermite .and. self%bc2 == sll_p_hermite) then
      allocate( val_corners (deg1/2, deg2/2, 4) )
      s1 = 1-modulo(deg1,2) ! shift = 1 for even order, 0 for odd order
      s2 = 1-modulo(deg2,2) ! shift = 1 for even order, 0 for odd order
      do j1 = 1, deg1/2
        do j2 = 1, deg2/2
          val_corners(j1,j2,1) = profile_2d % eval( info%x1_min, info%x2_min, diff_x1=j1-s1, diff_x2=j2-s2 )
          val_corners(j1,j2,2) = profile_2d % eval( info%x1_max, info%x2_min, diff_x1=j1-s1, diff_x2=j2-s2 )
          val_corners(j1,j2,3) = profile_2d % eval( info%x1_min, info%x2_max, diff_x1=j1-s1, diff_x2=j2-s2 )
          val_corners(j1,j2,4) = profile_2d % eval( info%x1_max, info%x2_max, diff_x1=j1-s1, diff_x2=j2-s2 )
        end do
      end do
    end if

    end associate ! tau1, tau2

    ! Compute 2D spline that interpolates analytical 2D profile at points above
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    call sll_s_set_time_mark( t2 )

    ! Hermite - other
    if (self%bc1 == sll_p_hermite .and. self%bc2 /= sll_p_hermite) then

      call sll_s_bspline_2d_compute_interpolant( self % bspline_2d, self % gtau, &
        val1_min = val1_min, &
        val1_max = val1_max )

    ! other - Hermite
    else if (self%bc1 /= sll_p_hermite .and. self%bc2 == sll_p_hermite) then

      call sll_s_bspline_2d_compute_interpolant( self % bspline_2d, self % gtau, &
        val2_min = val2_min, &
        val2_max = val2_max )

    ! Hermite - Hermite
    else if (self%bc1 == sll_p_hermite .and. self%bc2 == sll_p_hermite) then

      call sll_s_bspline_2d_compute_interpolant( self % bspline_2d, self % gtau, &
        val1_min    = val1_min, &
        val1_max    = val1_max, &
        val2_min    = val2_min, &
        val2_max    = val2_max, &
        val_corners = val_corners )

    ! other - other
    else

      call sll_s_bspline_2d_compute_interpolant( self % bspline_2d, self % gtau )

    end if

    call sll_s_set_time_mark( t3 )

    ! Deallocate local arrays
    ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if (allocated(val1_min))    deallocate( val1_min )
    if (allocated(val1_max))    deallocate( val1_max )
    if (allocated(val2_min))    deallocate( val2_min )
    if (allocated(val2_max))    deallocate( val2_max )
    if (allocated(val_corners)) deallocate( val_corners )

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

    class(t_bspline_2d_test_facility), intent(inout) :: self
    logical                          , intent(  out) :: equiv(3)

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
    !   . sll_f_bspline_2d_eval
    !   . sll_s_bspline_2d_eval_array
    call sll_s_set_time_mark( t0 )
    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_bspline_2d_eval( self%bspline_2d, x1(i1,i2), x2(i1,i2) )
      end do
    end do
    call sll_s_set_time_mark( t1 )
    call sll_s_bspline_2d_eval_array( self%bspline_2d, x1, x2, ya )
    call sll_s_set_time_mark( t2 )

    self%time_eval       = sll_f_time_elapsed_between( t0, t1 ) / npts
    self%time_eval_array = sll_f_time_elapsed_between( t1, t2 ) / npts
    equiv(1) = all( y == ya )

    ! Compare:
    !   . sll_f_bspline_2d_eval_deriv_x1
    !   . sll_s_bspline_2d_eval_array_deriv_x1
    call sll_s_set_time_mark( t0 )
    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_bspline_2d_eval_deriv_x1( self%bspline_2d, x1(i1,i2), x2(i1,i2) )
      end do
    end do
    call sll_s_set_time_mark( t1 )
    call sll_s_bspline_2d_eval_array_deriv_x1( self%bspline_2d, x1, x2, ya )
    call sll_s_set_time_mark( t2 )

    self%time_eval_diff1       = sll_f_time_elapsed_between( t0, t1 ) / npts
    self%time_eval_diff1_array = sll_f_time_elapsed_between( t1, t2 ) / npts
    equiv(2) = all( y == ya )

    ! Compare:
    !   . sll_f_bspline_2d_eval_deriv_x2
    !   . sll_s_bspline_2d_eval_array_deriv_x2
    call sll_s_set_time_mark( t0 )
    do i2 = 1, n2
      do i1 = 1, n1
        y(i1,i2) = sll_f_bspline_2d_eval_deriv_x2( self%bspline_2d, x1(i1,i2), x2(i1,i2) )
      end do
    end do
    call sll_s_set_time_mark( t1 )
    call sll_s_bspline_2d_eval_array_deriv_x2( self%bspline_2d, x1, x2, ya )
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

    class(t_bspline_2d_test_facility), intent(in   ) :: self
    real(wp)                         , intent(  out) :: max_norm_error

    integer  :: i1, i2
    real(wp) :: error

    ! Get spline interpolation points
    associate( tau1 => self % bspline_2d % bs1 % tau, &
               tau2 => self % bspline_2d % bs2 % tau )

      ! Evaluate 2D spline at interpolation points:
      ! interpolation values should be obtained
      max_norm_error = 0.0_wp
      do i2 = 1, size( tau2 )
        do i1 = 1, size( tau1 )
          error = self % gtau(i1,i2) &
                - sll_f_bspline_2d_eval( self % bspline_2d, tau1(i1), tau2(i2) )
          max_norm_error = max( max_norm_error, abs( error ) )
        end do
      end do

    end associate ! tau1, tau2

  end subroutine evaluate_at_interpolation_points

  !-----------------------------------------------------------------------------
  subroutine evaluate_on_2d_grid( self, x1, x2, max_norm_error )

    class(t_bspline_2d_test_facility), intent(inout) :: self
    real(wp)                         , intent(in   ) :: x1(:,:)
    real(wp)                         , intent(in   ) :: x2(:,:)
    real(wp)                         , intent(  out) :: max_norm_error

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
    call sll_s_bspline_2d_eval_array( self%bspline_2d, x1, x2, y )

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

    class(t_bspline_2d_test_facility), intent(inout) :: self
    real(wp)                         , intent(in   ) :: x1(:,:)
    real(wp)                         , intent(in   ) :: x2(:,:)
    real(wp)                         , intent(  out) :: max_norm_error_diff_x1
    real(wp)                         , intent(  out) :: max_norm_error_diff_x2

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

    call sll_s_bspline_2d_eval_array_deriv_x1( self%bspline_2d, x1, x2, y )

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

    call sll_s_bspline_2d_eval_array_deriv_x2( self%bspline_2d, x1, x2, y )

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

    class(t_bspline_2d_test_facility), intent(inout) :: self

    ! Free spline memory
    call sll_s_bspline_2d_free( self % bspline_2d )

  end subroutine free

end module m_test_bsplines_2d
