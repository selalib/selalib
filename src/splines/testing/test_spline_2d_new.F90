program test_spline_2d_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

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

  use m_test_spline_2d_new, only: &
    t_spline_2d_test_facility

  use m_analytical_profiles_2d, only: &
    t_profile_2d_info, &
    c_analytical_profile_2d, &
    t_analytical_profile_2d_cos_cos, &
    t_analytical_profile_2d_poly

  use m_splines_error_bounds, only: &
    sll_f_spline_2d_error_bound, &
    sll_f_spline_2d_error_bounds_on_grad

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  ! Local variables
  type(t_profile_2d_info)               :: pinfo
  type(t_spline_2d_test_facility)       :: test_facility
  type(t_analytical_profile_2d_cos_cos) :: profile_2d_cos_cos
  type(t_analytical_profile_2d_poly)    :: profile_2d_poly
  integer , allocatable                 :: bc_kinds(:)
  integer , allocatable                 :: nx_list(:)
  real(wp), allocatable                 :: breaks1(:)
  real(wp), allocatable                 :: breaks2(:)

  ! Evaluation grid
  real(wp), allocatable ::  grid_x1(:)
  real(wp), allocatable ::  grid_x2(:)
  real(wp), allocatable :: mgrid_x1(:,:)
  real(wp), allocatable :: mgrid_x2(:,:)

  ! Parameters for uniform / non-uniform grid
  logical  :: uniform
  real(wp) :: grid_perturbation

  integer  :: deg1
  integer  :: deg2
  integer  :: ncells (2)
  integer  :: bc_xmin(2)
  integer  :: bc_xmax(2)
  integer  :: i1, i2
  integer  :: j1, j2
  integer  :: grid_dim(2)
  real(wp) :: tol
  real(wp) :: tols_grad(2)
  real(wp) :: tol_diff_x1
  real(wp) :: tol_diff_x2
  logical  :: passed(3)
  logical  :: success
  logical  :: success_diff_x1
  logical  :: success_diff_x2
  real(wp) :: max_norm_error        , max_norm_profile
  real(wp) :: max_norm_error_diff_x1, max_norm_profile_diff_x1
  real(wp) :: max_norm_error_diff_x2, max_norm_profile_diff_x2
  real(wp) :: dx1
  real(wp) :: dx2
  integer  :: k
  integer  :: cos_n1, cos_n2
  real(wp) :: cos_c1, cos_c2

  ! Read from standard input
  call process_args( uniform, grid_perturbation )

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 1: Evaluate spline at interpolation points (error should be zero)
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  write(*,*)
  write(*,'(a)') '------------------------------------------------------------------------'
  write(*,'(a)') ' TEST 1: evaluate spline at interpolation points (error should be zero) '
  write(*,'(a)') '------------------------------------------------------------------------'
  write(*,*)

  ! Test tolerance: very small!
  tol = 1e-14_wp

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(3) )
  bc_kinds(:) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  ! Choose number of cells in tensor grid
  ncells = [11, 36]

  ! Print constant parameters
  write(*,'(a)') 'Profile: f(x1,x2) = cos( 2*pi*x1 ) * cos( 2*pi*x2 )'
  write(*,*)
  write(*,'(a,es10.1)') 'Relative error tolerance: tol = ', tol
  write(*,'(a,i4)')     'Number of cells in grid : nx1 = ', ncells(1)
  write(*,'(a,i4)')     '                          nx2 = ', ncells(2)
  write(*,*)
  write(*,'(a)') "Input:"
  write(*,'(a)') '  . deg1 = spline degree in x1 direction'
  write(*,'(a)') '  . deg2 = spline degree in x2 direction'
  write(*,'(a)') '  .  bc1 = boundary conditions in x1 direction [P|H|G]'
  write(*,'(a)') '  .  bc2 = boundary conditions in x2 direction [P|H|G]'
  write(*,*)
  write(*,'(a)') "Output:"
  write(*,'(a)') '  .  error  = relative max-norm of error'
  write(*,'(a)') '  .  passed = "OK" if error <= tol, "FAIL" otherwise'
  write(*,*)
  write(*,'(a)') "Boundary conditions:"
  write(*,'(a)') '  .  P = periodic'
  write(*,'(a)') '  .  H = Hermite'
  write(*,'(a)') '  .  G = Greville'
  write(*,*)

  ! Print table header
  write(*, '(4a6,2a10)') "deg1", "deg2", "bc1", "bc2", "error", "passed"

  ! Initialize profile
  call profile_2d_cos_cos % init()

  ! Extract information about 2D analytical profile
  call profile_2d_cos_cos % get_info( pinfo )

  ! Estimate max-norm of profile (needed to compute relative error)
  max_norm_profile = profile_2d_cos_cos % max_norm()

  if (uniform) then
    allocate( breaks1(0) )
    allocate( breaks2(0) )
  else
    allocate( breaks1(ncells(1)+1) )
    allocate( breaks2(ncells(2)+1) )
    call generate_non_uniform_breaks( pinfo%x1_min, pinfo%x1_max, ncells(1), grid_perturbation, breaks1 )
    call generate_non_uniform_breaks( pinfo%x2_min, pinfo%x2_max, ncells(2), grid_perturbation, breaks2 )
  end if

  ! Initialize 'PASSED/FAILED' condition
  passed(1) = .true.

  ! Cycle over spline degree
  do deg1 = 1, 9
    do deg2 = 1, 9

      ! Cycle over all kinds of boundary conditions at x1_min
      do i1 = 1, size( bc_kinds )
        bc_xmin(1) = bc_kinds(i1)
        if (bc_xmin(1) == sll_p_periodic .and. (.not. pinfo % x1_periodic)) cycle

        ! Cycle over all kinds of boundary conditions at x1_max
        do j1 = 1, size( bc_kinds )
          bc_xmax(1) = bc_kinds(j1)
          if ( any( [bc_xmin(1),bc_xmax(1)] == sll_p_periodic ) &
              .and. bc_xmin(1) /= bc_xmax(1) ) cycle

          ! Cycle over all kinds of boundary conditions in x2_min
          do i2 = 1, size( bc_kinds )
            bc_xmin(2) = bc_kinds(i2)
            if (bc_xmin(2) == sll_p_periodic .and. (.not. pinfo % x2_periodic)) cycle

            ! Cycle over all kinds of boundary conditions in x2_max
            do j2 = 1, size( bc_kinds )
              bc_xmax(2) = bc_kinds(j2)
              if ( any( [bc_xmin(2),bc_xmax(2)] == sll_p_periodic ) &
                  .and. bc_xmin(2) /= bc_xmax(2) ) cycle

              ! Initialize test facility
              call test_facility % init( &
                profile_2d = profile_2d_cos_cos, &
                degree     = [deg1, deg2], &
                ncells     = ncells (1:2), &
                bc_xmin    = bc_xmin(1:2), &
                bc_xmax    = bc_xmax(1:2), &
                breaks1    = breaks1(:)  , &
                breaks2    = breaks2(:)  )

              ! Run tests
              call test_facility % evaluate_at_interpolation_points( max_norm_error )

              ! Calculate relative error norm from absolute one, check tolerance
              max_norm_error = max_norm_error / max_norm_profile
              success = (max_norm_error <= tol)

              ! Print test report to terminal on a single line
              write(*,'(2i6)', advance='no') deg1, deg2
              write(*,'(2a6)', advance='no') &
                bc_to_char( bc_xmin(1) ) // "-" // bc_to_char( bc_xmax(1) ), &
                bc_to_char( bc_xmin(2) ) // "-" // bc_to_char( bc_xmax(2) )
              write(*,'(es10.1,a8)') max_norm_error, trim( success_to_string( success ))

              ! Free memory
              call test_facility % free()

              ! Update 'PASSED/FAILED' condition
              passed(1) = (passed(1) .and. success)

            end do  ! bc_xmax(2)
          end do  ! bc_xmin(2)
        end do  ! bc_xmax(1)
      end do  ! bc_xmin(1)
    end do  ! deg2
  end do  ! deg1

  ! Deallocate local arrays
  deallocate( bc_kinds )
  deallocate( breaks1  )
  deallocate( breaks2  )

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 2: Spline should represent polynomial profiles exactly
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  write(*,*)
  write(*,'(a)') '---------------------------------------------------------------------------'
  write(*,'(a)') ' TEST 2: spline of order (m,n) exactly represents polynomial of same order '
  write(*,'(a)') '---------------------------------------------------------------------------'
  write(*,*)

  ! Test tolerance: very small!
  tol         = 1e-14_wp
  tol_diff_x1 = 2e-13_wp ! larger for derivatives
  tol_diff_x2 = 2e-13_wp ! larger for derivatives

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(2) )
  bc_kinds(:) = [sll_p_hermite, sll_p_greville]

  ! Choose number of cells in tensor grid
  ncells = [22, 17]

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = [20, 20]

  ! Print constant parameters
  write(*,'(a,es10.1)') 'Relative error tolerance for f     : tol         = ', tol
  write(*,'(a,es10.1)') 'Relative error tolerance for ∂f/∂x1: tol_diff_x1 = ', tol_diff_x1
  write(*,'(a,es10.1)') 'Relative error tolerance for ∂f/∂x2: tol_diff_x2 = ', tol_diff_x2
  write(*,'(a,i4)')     'Number of cells in grid: nx1 = ', ncells(1)
  write(*,'(a,i4)')     '                         nx2 = ', ncells(2)
  write(*,'(a,i4,a,i4,a)')'Number of evaluation points (uniform grid): [', grid_dim(1), ',', grid_dim(2),']'
  write(*,*)
  write(*,'(a)') "Input:"
  write(*,'(a)') '  . deg1 = spline degree in x1 direction'
  write(*,'(a)') '  . deg2 = spline degree in x2 direction'
  write(*,'(a)') '  .  bc1 = boundary conditions in x1 direction [H|G]'
  write(*,'(a)') '  .  bc2 = boundary conditions in x2 direction [H|G]'
  write(*,*)
  write(*,'(a)') "Output:"
  write(*,'(a)') '  .  err     = relative max-norm of error on f'
  write(*,'(a)') '  .  err_dx1 = relative max-norm of error on ∂f/∂x1'
  write(*,'(a)') '  .  err_dx2 = relative max-norm of error on ∂f/∂x2'
  write(*,'(a)') '  .  passed  = "OK" if all errors <= tol, "FAIL" otherwise'
  write(*,*)
  write(*,'(a)') "Boundary conditions:"
  write(*,'(a)') '  .  H = Hermite'
  write(*,'(a)') '  .  G = Greville'
  write(*,*)

  ! Print table header
  write(*, '(4a6,3a10,a10)') "deg1", "deg2", "bc1", "bc2", &
    "err", "err_dx1", "err_dx2", "passed"

  call profile_2d_poly % init( 0, 0 )  ! dummy, only need domain size
  call profile_2d_poly % get_info( pinfo )

  ! Create 1D uniform grids of evaluation points
  allocate( grid_x1 (grid_dim(1)) )
  allocate( grid_x2 (grid_dim(2)) )
  call sll_s_new_array_linspace( grid_x1, pinfo%x1_min, pinfo%x1_max, step=dx1 )
  call sll_s_new_array_linspace( grid_x2, pinfo%x2_min, pinfo%x2_max, step=dx2 )

  ! Create 2D meshgrids from 1D grids (this could be done with a subroutine)
  allocate( mgrid_x1 (grid_dim(1),grid_dim(2)) )
  allocate( mgrid_x2 (grid_dim(1),grid_dim(2)) )
  do i2 = 1, grid_dim(2)
    do i1 = 1, grid_dim(1)
      mgrid_x1(i1,i2) = grid_x1(i1)
      mgrid_x2(i1,i2) = grid_x2(i2)
    end do
  end do

  if (uniform) then
    allocate( breaks1(0) )
    allocate( breaks2(0) )
  else
    allocate( breaks1(ncells(1)+1) )
    allocate( breaks2(ncells(2)+1) )
    call generate_non_uniform_breaks( pinfo%x1_min, pinfo%x1_max, ncells(1), grid_perturbation, breaks1 )
    call generate_non_uniform_breaks( pinfo%x2_min, pinfo%x2_max, ncells(2), grid_perturbation, breaks2 )
  end if

  ! Initialize 'PASSED/FAILED' condition
  passed(2) = .true.

  do deg1 = 1, 9
    do deg2 = 1, 9

      ! Initialize polynomial profile with degree (deg1,deg2)
      call profile_2d_poly % init( deg1, deg2 )

      ! Estimate max-norm of profile (needed to compute relative error)
      max_norm_profile         = profile_2d_poly % max_norm()
      max_norm_profile_diff_x1 = profile_2d_poly % max_norm( diff_x1=1 )
      max_norm_profile_diff_x2 = profile_2d_poly % max_norm( diff_x2=1 )

      ! Cycle over all kinds of boundary conditions at x1_min
      do i1 = 1, size( bc_kinds )
        bc_xmin(1) = bc_kinds(i1)

        ! Cycle over all kinds of boundary conditions at x1_max
        do j1 = 1, size( bc_kinds )
          bc_xmax(1) = bc_kinds(j1)

          ! Cycle over all kinds of boundary conditions in x2_min
          do i2 = 1, size( bc_kinds )
            bc_xmin(2) = bc_kinds(i2)

            ! Cycle over all kinds of boundary conditions in x2_max
            do j2 = 1, size( bc_kinds )
              bc_xmax(2) = bc_kinds(j2)

              ! Initialize test facility
              call test_facility % init( &
                profile_2d = profile_2d_poly, &
                degree     = [deg1, deg2], &
                ncells     = ncells (1:2), &
                bc_xmin    = bc_xmin(1:2), &
                bc_xmax    = bc_xmax(1:2), &
                breaks1    = breaks1(:)  , &
                breaks2    = breaks2(:) )

              ! Run tests
              call test_facility % evaluate_on_2d_grid( mgrid_x1, mgrid_x2, max_norm_error )
              call test_facility % evaluate_grad_on_2d_grid( mgrid_x1, mgrid_x2, &
                max_norm_error_diff_x1, max_norm_error_diff_x2 )

              ! Calculate relative error norms from absolute ones
              max_norm_error         = max_norm_error         / max_norm_profile
              max_norm_error_diff_x1 = max_norm_error_diff_x1 / max_norm_profile_diff_x1
              max_norm_error_diff_x2 = max_norm_error_diff_x2 / max_norm_profile_diff_x2

              ! Check tolerances
              success         = (max_norm_error         <= tol        )
              success_diff_x1 = (max_norm_error_diff_x1 <= tol_diff_x1)
              success_diff_x2 = (max_norm_error_diff_x2 <= tol_diff_x2)

              ! Keep single success condition
              success = success .and. success_diff_x1 .and. success_diff_x2

              ! Print test report to terminal on a single line
              write(*,'(2i6)'    , advance='no') deg1, deg2
              write(*,'(2a6)', advance='no') &
                bc_to_char( bc_xmin(1) ) // "-" // bc_to_char( bc_xmax(1) ), &
                bc_to_char( bc_xmin(2) ) // "-" // bc_to_char( bc_xmax(2) )
              write(*,'(3es10.1)', advance='no') &
                max_norm_error, max_norm_error_diff_x1, max_norm_error_diff_x2
              write(*,'(a8)') trim( success_to_string( success ))

              ! Free memory
              call test_facility % free()

              ! Update 'PASSED/FAILED' condition
              passed(2) = (passed(2) .and. success)

            end do  ! bc_xmax(2)
          end do  ! bc_xmin(2)
        end do  ! bc_xmax(1)
      end do  ! bc_xmin(1)
    end do  ! deg2
  end do  ! deg1

  ! Deallocate local arrays
  deallocate( bc_kinds )
  deallocate(  grid_x1 )
  deallocate(  grid_x2 )
  deallocate( mgrid_x1 )
  deallocate( mgrid_x2 )
  deallocate( breaks1  )
  deallocate( breaks2  )

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 3: convergence analysis on cos*cos profile (with absolute error bound)
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  write(*,*)
  write(*,'(a)') '-----------------------------------------------------------------------------'
  write(*,'(a)') ' TEST 3: convergence analysis on cos*cos profile (with absolute error bound) '
  write(*,'(a)') '-----------------------------------------------------------------------------'
  write(*,*)

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(3) )
  bc_kinds(:) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  allocate( nx_list(5) )
  nx_list(:) = [10, 20, 40, 80, 160]

  ! Choose parameters of cos*cos profile to be interpolated
  cos_n1 = 3
  cos_n2 = 3
  cos_c1 = 0.3_wp
  cos_c2 = 0.7_wp

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = [20,20]

  ! Print constant parameters
  write(*,'(a)') 'Profile: f(x1,x2) = cos( 2*pi * (n1*x1+c1) ) * cos( 2*pi * (n2*x2+c2) )'
  write(*,'(a,i2)'    ) '  . n1 = ', cos_n1
  write(*,'(a,i2)'    ) '  . n2 = ', cos_n2
  write(*,'(a,f17.15)') '  . c1 = ', cos_c1
  write(*,'(a,f17.15)') '  . c2 = ', cos_c2
  write(*,*)
  write(*,'(a,i4,a,i4,a)')'Number of evaluation points (uniform grid): [', grid_dim(1), ',', grid_dim(2),']'
  write(*,*)
  write(*,'(a)') "Input:"
  write(*,'(a)') '  . deg = spline degree (same along x1 and x2 directions)'
  write(*,'(a)') '  . bc1 = boundary conditions in x1 direction [P|H|G]'
  write(*,'(a)') '  . bc2 = boundary conditions in x2 direction [P|H|G]'
  write(*,'(a)') '  . nx  = number of cells along each direction (nx1=nx2)'
  write(*,*)
  write(*,'(a)') "Output:"
  write(*,'(a)') '  . tol     = tolerance for f'
  write(*,'(a)') '  . err     = max-norm of error on f'
  write(*,'(a)') '  . tol_dx1 = tolerance for ∂f/∂x1'
  write(*,'(a)') '  . err_dx1 = max-norm of error on ∂f/∂x1'
  write(*,'(a)') '  . tol_dx2 = tolerance for ∂f/∂x2'
  write(*,'(a)') '  . err_dx2 = max-norm of error on ∂f/∂x2'
  write(*,'(a)') '  . passed  = "OK" if all errors <= tol, "FAIL" otherwise'
  write(*,*)
  write(*,'(a)') "Timing [s]:"
  write(*,'(a)') '  . t_init  = initialization of spline object'
  write(*,'(a)') '  . t_comp  = calculation of spline coefficients for interpolation'
  write(*,'(a)') '  . t_eval  = evaluation of function value S(x1,x2)'
  write(*,'(a)') '  . t_dx1   = evaluation of x1-derivative ∂S(x1,x2)/∂x1'
  write(*,'(a)') '  . t_dx2   = evaluation of x2-derivative ∂S(x1,x2)/∂x2'
  write(*,*)
  write(*,'(a)') "Boundary conditions:"
  write(*,'(a)') '  . P = periodic'
  write(*,'(a)') '  . H = Hermite'
  write(*,'(a)') '  . G = Greville'
  write(*,*)

  ! Print table header
  write(*, '(a3,3a6,*(a10))') &
    "deg", "bc1", "bc2", "nx", &
    "tol", "err", "tol_dx1", "err_dx1", "tol_dx2", "err_dx2", "passed", &
    "t_init", "t_comp", "t_eval", "t_dx1", "t_dx2"

  ! Initialize profile
  call profile_2d_cos_cos % init( cos_n1, cos_n2, cos_c1, cos_c2 )

  ! Extract information about 2D analytical profile
  call profile_2d_cos_cos % get_info( pinfo )

  ! Estimate max-norm of profile (needed to compute relative error)
  max_norm_profile         = profile_2d_cos_cos % max_norm()
  max_norm_profile_diff_x1 = profile_2d_cos_cos % max_norm( diff_x1=1 )
  max_norm_profile_diff_x2 = profile_2d_cos_cos % max_norm( diff_x2=1 )

  ! Create 1D uniform grids of evaluation points
  allocate( grid_x1 (grid_dim(1)) )
  allocate( grid_x2 (grid_dim(2)) )
  call sll_s_new_array_linspace( grid_x1, pinfo%x1_min, pinfo%x1_max, step=dx1 )
  call sll_s_new_array_linspace( grid_x2, pinfo%x2_min, pinfo%x2_max, step=dx2 )

  ! Create 2D meshgrids from 1D grids (this could be done with a subroutine)
  allocate( mgrid_x1 (grid_dim(1),grid_dim(2)) )
  allocate( mgrid_x2 (grid_dim(1),grid_dim(2)) )
  do i2 = 1, grid_dim(2)
    do i1 = 1, grid_dim(1)
      mgrid_x1(i1,i2) = grid_x1(i1)
      mgrid_x2(i1,i2) = grid_x2(i2)
    end do
  end do

  ! Initialize 'PASSED/FAILED' condition
  passed(3) = .true.

  do deg1 = 1, 9

    ! Use same spline degree along both directions
    deg2 = deg1

    ! Cycle over all kinds of boundary conditions at x1_min
    do i1 = 1, size( bc_kinds )
      bc_xmin(1) = bc_kinds(i1)

      ! Cycle over all kinds of boundary conditions at x1_max
      do j1 = 1, size( bc_kinds )
        bc_xmax(1) = bc_kinds(j1)
        if ( any( [bc_xmin(1),bc_xmax(1)] == sll_p_periodic ) &
            .and. bc_xmin(1) /= bc_xmax(1) ) cycle

        ! Cycle over all kinds of boundary conditions in x2_min
        do i2 = 1, size( bc_kinds )
          bc_xmin(2) = bc_kinds(i2)

          ! Cycle over all kinds of boundary conditions in x2_max
          do j2 = 1, size( bc_kinds )
            bc_xmax(2) = bc_kinds(j2)
            if ( any( [bc_xmin(2),bc_xmax(2)] == sll_p_periodic ) &
                .and. bc_xmin(2) /= bc_xmax(2) ) cycle

            ! Convergence analysis: increase number of knots
            do k = 1, size( nx_list )
              ncells(1) = nx_list(k)

              ! Use same number of knots along both directions
              ncells(2) = ncells(1)

              if (uniform) then
                allocate( breaks1(0) )
                allocate( breaks2(0) )
                dx1 = (pinfo%x1_max - pinfo%x1_min) / ncells(1)
                dx2 = (pinfo%x2_max - pinfo%x2_min) / ncells(2)
              else
                allocate( breaks1(ncells(1)+1) )
                allocate( breaks2(ncells(2)+1) )
                call generate_non_uniform_breaks( pinfo%x1_min, pinfo%x1_max, ncells(1), grid_perturbation, breaks1 )
                call generate_non_uniform_breaks( pinfo%x2_min, pinfo%x2_max, ncells(2), grid_perturbation, breaks2 )
                dx1 = maxval( breaks1(2:ncells(1)+1)-breaks1(1:ncells(1)) )
                dx2 = maxval( breaks2(2:ncells(2)+1)-breaks2(1:ncells(2)) )
              end if

              ! Initialize test facility
              call test_facility % init( &
                profile_2d = profile_2d_cos_cos, &
                degree     = [deg1, deg2], &
                ncells     = ncells (1:2), &
                bc_xmin    = bc_xmin(1:2), &
                bc_xmax    = bc_xmax(1:2), &
                breaks1    = breaks1(:)  , &
                breaks2    = breaks2(:)  )

              deallocate( breaks1 )
              deallocate( breaks2 )

              ! Determine error tolerances
              tol       = sll_f_spline_2d_error_bound         ( profile_2d_cos_cos, dx1, dx2, deg1, deg2 )
              tols_grad = sll_f_spline_2d_error_bounds_on_grad( profile_2d_cos_cos, dx1, dx2, deg1, deg2 )

              ! Increase absolute tolerances if below machine precision
              ! NOTE: derivatives more sensitive to roundoff
              tol         = max( 1e-14_wp*max_norm_profile        , tol          )
              tol_diff_x1 = max( 1e-12_wp*max_norm_profile_diff_x1, tols_grad(1) )
              tol_diff_x2 = max( 1e-12_wp*max_norm_profile_diff_x2, tols_grad(2) )

              ! Run tests
              ! TODO: print numerical order of accuracy
              call test_facility % evaluate_on_2d_grid( mgrid_x1, mgrid_x2, max_norm_error )
              call test_facility % evaluate_grad_on_2d_grid( mgrid_x1, mgrid_x2, &
                max_norm_error_diff_x1, max_norm_error_diff_x2 )

              ! Check tolerances
              success         = (max_norm_error         <= tol        )
              success_diff_x1 = (max_norm_error_diff_x1 <= tol_diff_x1)
              success_diff_x2 = (max_norm_error_diff_x2 <= tol_diff_x2)

              ! Keep single success condition
              success = success .and. success_diff_x1 .and. success_diff_x2

              ! Print test report to terminal on a single line
              write(*,'(i3)'     , advance='no') deg1
              write(*,'(2a6)'    , advance='no') &
                bc_to_char( bc_xmin(1) ) // "-" // bc_to_char( bc_xmax(1) ), &
                bc_to_char( bc_xmin(2) ) // "-" // bc_to_char( bc_xmax(2) )
              write(*,'(2i6)'    , advance='no') ncells(1)
              write(*,'(2es10.1)', advance='no') tol        , max_norm_error
              write(*,'(2es10.1)', advance='no') tol_diff_x1, max_norm_error_diff_x1
              write(*,'(2es10.1)', advance='no') tol_diff_x2, max_norm_error_diff_x2
              write(*,'(a10)'    , advance='no') success_to_string( success )
              write(*,'(5es10.1)') &
                test_facility % time_init               , &
                test_facility % time_compute_interpolant, &
                test_facility % time_eval_array         , &
                test_facility % time_eval_diff1_array   , &
                test_facility % time_eval_diff2_array

              ! Free memory
              call test_facility % free()

              ! Update 'PASSED/FAILED' condition
              passed(3) = (passed(3) .and. success)

            end do  ! nx1
            write(*,*)

          end do  ! bc_xmax(2)
        end do  ! bc_xmin(2)
      end do  ! bc_xmax(1)
    end do  ! bc_xmin(1)
  end do   ! deg1=deg2

  ! Deallocate local arrays
  deallocate( bc_kinds  )
  deallocate( nx_list   )
  deallocate(  grid_x1  )
  deallocate(  grid_x2  )
  deallocate( mgrid_x1  )
  deallocate( mgrid_x2  )

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  write(*,'(a)') '----------------'
  write(*,'(a)') 'Test 1: '// success_to_string( passed(1) )
  write(*,'(a)') 'Test 2: '// success_to_string( passed(2) )
  write(*,'(a)') 'Test 3: '// success_to_string( passed(3) )
  write(*,'(a)') '----------------'
  write(*,*)

  ! CTest key
  if (all( passed )) then
    write(*,'(a)') "CTEST: PASSED"
  else
    write(*,'(a)') "CTEST: FAILED"
  end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  pure function bc_to_char( bc ) result( c )
    integer, intent(in) :: bc
    character(len=1) :: c

    select case (bc)
    case (sll_p_periodic)
      c = 'P'
    case (sll_p_hermite)
      c = 'H'
    case (sll_p_greville)
      c = 'G'
    case default
      c = '?'
    end select

  end function bc_to_char

  !-----------------------------------------------------------------------------
  pure function success_to_string( success ) result( str )
    logical, intent(in) :: success
    character(len=4) :: str

    if (success) then
      str = 'OK'
    else
      str = 'FAIL'
    end if

  end function success_to_string

  !-----------------------------------------------------------------------------
  pure function factorial( n ) result( f )
    integer, intent(in) :: n
    real(wp) :: f

    f = gamma( real( n+1, wp ) )

  end function factorial

  !-----------------------------------------------------------------------------
  subroutine process_args( uniform, grid_perturbation )
    logical , intent(out) :: uniform
    real(wp), intent(out) :: grid_perturbation

    integer :: argc
    character(len=32) :: val
    integer :: length
    integer :: status

    argc = command_argument_count()

    select case (argc)

    case(0)
      SLL_ERROR("process_args","command line arguments: -u | -n [ real number in [0,1) ]")

    case(1)
      call get_command_argument(1,val,length,status)
      if (trim(val) == "-u") then
        uniform = .true.
      else if (trim(val) == "-n") then
        uniform = .false.
        grid_perturbation = 0.0_wp
      else
        SLL_ERROR("process_args","command line arguments: -u | -n [ real number in [0,1) ]")
      end if

    case(2)
      call get_command_argument(1,val,length,status)
      if (trim(val) == "-u") then
        SLL_ERROR("process_args","command line arguments: -u | -n [ real number in [0,1) ]")
      else if (trim(val) == "-n") then
        uniform = .false.
      else
        SLL_ERROR("process_args","command line arguments: -u | -n [ real number in [0,1) ]")
      end if
      call get_command_argument(2,val,length,status)
      read(val,*) grid_perturbation

    case default
      SLL_ERROR("process_args","command line arguments: -u | -n [grid perturbation]")

    end select

  end subroutine process_args

  !-----------------------------------------------------------------------------
  subroutine generate_non_uniform_breaks( xmin, xmax, ncells, grid_perturbation, breaks )
    real(wp), intent(in   ) :: xmin
    real(wp), intent(in   ) :: xmax
    integer , intent(in   ) :: ncells
    real(wp), intent(in   ) :: grid_perturbation
    real(wp), intent(  out) :: breaks(:)

    integer  :: i
    real(wp) :: r
    real(wp) :: a, b

    ! Generate breakpoints by applying random noise onto regular grid
    associate( dx => (xmax-xmin)/ncells )
      do i = 1, ncells+1
        call random_number( r ) !  0.0 <= r < 1.0
        r = r - 0.5_wp          ! -0.5 <= r < 0.5
        breaks(i) = xmin + ( real(i-1,wp) + grid_perturbation * r ) * dx
      end do
    end associate

    ! Linearly transform coordinates in order to match nominal xmin and xmax
    associate( y    => breaks(:), &
               ymin => breaks(1), &
               ymax => breaks(ncells+1) )

      a = (xmin-xmax)/(ymin-ymax)
      b = (xmax*ymin-xmin*ymax)/(ymin-ymax)

      breaks(1) = xmin
      do i = 2, ncells
        breaks(i) =  a * y(i) + b
      end do
      breaks(ncells+1) = xmax

    end associate

  end subroutine generate_non_uniform_breaks

end program test_spline_2d_new
