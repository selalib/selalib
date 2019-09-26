program test_polar_spline_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: &
    f64

  use sll_m_constants, only: &
    sll_p_twopi

  use sll_m_utilities, only: &
    sll_s_new_array_linspace

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite, &
    sll_p_greville

  use sll_m_bsplines, only: &
    sll_c_bsplines, &
    sll_s_bsplines_new, &
    sll_s_bsplines_new_2d_polar, &
    sll_s_bsplines_new_mirror_copy

  use sll_m_spline_2d, only: &
    sll_t_spline_2d

  use sll_m_polar_spline_interpolator_2d, only: &
    sll_t_polar_spline_interpolator_2d, &
    sll_s_polar_spline_2d_compute_num_cells

  use sll_m_timer, only: &
    sll_t_time_mark, &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_between

  use m_analytical_profiles_2d_base, only: &
    t_profile_2d_info, &
    c_analytical_profile_2d

  use m_analytical_profiles_2d_cos_cos, only: &
    t_analytical_profile_2d_cos_cos

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  class(sll_c_bsplines), allocatable       :: bsplines_x1
  class(sll_c_bsplines), allocatable       :: bsplines_x2
  type(sll_t_spline_2d)                    :: spline_2d
  type(sll_t_polar_spline_interpolator_2d) :: spline_interpolator_2d

  real(wp), allocatable ::  tau1(:)   ! interp. points, x1 coord
  real(wp), allocatable ::  tau2(:)   ! interp. points, x2 coord
  real(wp), allocatable :: gtau (:,:) ! profile values at tau

  real(wp) :: r_max
  real(wp) :: th_period
  integer  :: nipts (2)
  integer  :: degree(2)
  integer  :: bc_rmax
  integer  :: ncells(2)
  integer  :: n1, n2
  real(wp) :: c1, c2
  integer  :: i1, i2
  real(wp) :: x(2), eta(2)
  real(wp) :: error
  real(wp) :: max_error
  integer  :: ntest_pts(2)
  real(wp), allocatable :: test_grid_eta1(:)
  real(wp), allocatable :: test_grid_eta2(:)

  type(t_analytical_profile_2d_cos_cos) :: profile
  !-----------------------------------------------------------------------------

  ! Domain info
  r_max = 1.0_wp
  th_period = sll_p_twopi

  ! Spline info
  nipts   = [10, 20] * 8
  degree  = [ 7,  7]
  bc_rmax = sll_p_greville

  ! Profile info (cos-cos)
  n1 = 2
  n2 = 3
  c1 = 0.0_wp
  c2 = 0.0_wp

  ! Test grid
  ntest_pts = [11111, 555]

  !-----------------------------------------------------------------------------
  call sll_s_polar_spline_2d_compute_num_cells( degree, bc_rmax, nipts, ncells )

  write(*,*) ""
  write(*,*) ncells

  ! Create uniform B-splines along x1 and x2 on polar grid
  call sll_s_bsplines_new_2d_polar( &
    bsplines_radial  = bsplines_x1, &
    bsplines_angular = bsplines_x2, &
    degree           = degree     , &
    ncells           = ncells     , &
    max_radius       = r_max      , &
    theta_lims       = [0.0_wp, th_period] )

  ! Initialize 2D polar spline
  call spline_2d % init( bsplines_x1, bsplines_x2 )

  ! Initialize 2D polar spline interpolator
  call spline_interpolator_2d % init( bsplines_x1, bsplines_x2, bc_rmax )

  ! Obtain interpolation points and allocate 2D array of values
  call spline_interpolator_2d % get_interp_points( tau1, tau2 )
  allocate( gtau (size(tau1),size(tau2)) )

  write(*,*) ""
  write(*,*) tau1(1:2)
!  write(*,*) ""
!  write(*,*) tau2

  ! Evaluate some analytical profile at the interpolation points
  call profile % init( n1, n2, c1, c2 )

  do i1 = 1, nipts(1)
    do i2 = 1, nipts(2)
      x = polar_to_cartesian( [tau1(i1), tau2(i2)] )
      gtau(i1,i2) = profile % eval( x(1), x(2) )
    end do
  end do

!  write(*,*) ""
!  write(*,*) gtau

  ! Compute spline that interpolates profile on grid
  call spline_interpolator_2d % compute_interpolant( spline_2d, gtau )

  write(*,*) ""
  write(*,*) "-----------------------------------------------------------------"
  write(*,*) " Evaluate spline at interpolation points "
  write(*,*) "-----------------------------------------------------------------"

  ! Evaluate spline at interpolation points
  max_error = 0.0_wp
  do i1 = 1, nipts(1)
    do i2 = 1, nipts(2)
      error = gtau(i1,i2) - spline_2d % eval( tau1(i1), tau2(i2) )
      max_error = max( max_error, abs( error ) )
    end do
  end do

  write(*,*) ""
  write(*,*) "Max error = ", max_error

  write(*,*) ""
  write(*,*) "-----------------------------------------------------------------"
  write(*,*) " Evaluate spline on test grid "
  write(*,*) "-----------------------------------------------------------------"

  ! Create test grid
  allocate( test_grid_eta1 (ntest_pts(1)) )
  allocate( test_grid_eta2 (ntest_pts(2)) )

  call sll_s_new_array_linspace( test_grid_eta1, 0.0_wp, r_max    , endpoint=.true.  )
  call sll_s_new_array_linspace( test_grid_eta2, 0.0_wp, th_period, endpoint=.false. )

  ! Evaluate spline on test grid
  max_error = 0.0_wp
  do i1 = 1, ntest_pts(1)
    do i2 = 1, ntest_pts(2)
      eta = [test_grid_eta1(i1), test_grid_eta2(i2)]
      x   = polar_to_cartesian( eta )
      error = profile % eval( x(1), x(2) ) - spline_2d % eval( eta(1), eta(2) )
      max_error = max( max_error, abs( error ) )
    end do
  end do

  write(*,*) ""
  write(*,*) "Max error = ", max_error
  write(*,*) ""

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  pure function polar_to_cartesian( eta ) result( x )
    real(wp), intent(in) :: eta(2)
    real(wp) :: x(2)

    associate( r => eta(1), theta => eta(2) )
      x(1) = r * cos( theta )
      x(2) = r * sin( theta )
    end associate

  end function polar_to_cartesian

end program test_polar_spline_2d
