program test_spline_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: &
    f64

  use sll_m_constants, only: &
    sll_p_twopi

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite, &
    sll_p_greville

  use m_test_spline_1d, only: &
    t_spline_1d_test_facility

  use m_analytical_profiles_1d, only: &
    t_profile_1d_info, &
    c_analytical_profile_1d, &
    t_analytical_profile_1d_cos, &
    t_analytical_profile_1d_poly

  use m_splines_error_bounds, only: &
    sll_f_spline_1d_error_bound, &
    sll_f_spline_1d_error_bound_on_deriv

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  ! Local variables
  type(t_profile_1d_info)               :: pinfo
  type(t_spline_1d_test_facility)       :: test_facility
  type(t_analytical_profile_1d_cos )    :: profile_1d_cos
  type(t_analytical_profile_1d_poly)    :: profile_1d_poly
  integer , allocatable                 :: bc_kinds(:)
  integer , allocatable                 :: nx_list(:)
  real(wp), allocatable                 :: grid(:)

  integer  :: degree
  integer  :: ncells
  integer  :: bc_xmin
  integer  :: bc_xmax
  integer  :: i, j, k

  integer  :: grid_dim
  real(wp) :: grid_dx

  integer  :: cos_n
  real(wp) :: cos_c
  real(wp) :: dx
  real(wp) :: tol
  real(wp) :: tol_diff
  real(wp) :: max_norm_profile
  real(wp) :: max_norm_profile_diff
  real(wp) :: max_norm_error
  real(wp) :: max_norm_error_diff
  logical  :: passed(1:3)
  logical  :: success
  logical  :: success_diff

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

  ! Choose number of cells in grid
  ncells = 12

  ! Print constant parameters
  write(*,'(a)') 'Profile: f(x) = cos( 2*pi*x )'
  write(*,*)
  write(*,'(a,es10.1)') 'Relative error tolerance: tol     = ', tol
  write(*,'(a,i4)')     'Number of cells in grid : ncells  = ', ncells
  write(*,*)
  write(*,'(a)') "Input:"
  write(*,'(a)') '  . deg   = spline degree'
  write(*,'(a)') '  . bcmin = boundary conditions at x=xmin [P|H|G]'
  write(*,'(a)') '  . bcmax = boundary conditions at x=xmax [P|H|G]'
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
  write(*, '(3a6,a10,a10)') "deg", "bcmin", "bcmax", "error", "passed"

  ! Initialize profile
  call profile_1d_cos % init()

  ! Extract information about 2D analytical profile
  call profile_1d_cos % get_info( pinfo )

  ! Estimate max-norm of profile (needed to compute relative error)
  max_norm_profile = profile_1d_cos % max_norm()

  ! Compute cell size in uniform grid
  dx = (pinfo%xmax-pinfo%xmin) / ncells

  ! Initialize 'PASSED/FAILED' condition
  passed(1) = .true.

  ! Cycle over spline degree
  do degree = 1, 9

    ! Cycle over all kinds of boundary conditions at x=xmin
    do i = 1, size( bc_kinds )
      bc_xmin = bc_kinds(i)
      if (bc_xmin == sll_p_periodic .and. (.not. pinfo % periodic)) cycle

      ! Cycle over all kinds of boundary conditions at x=xmax
      do j = 1, size( bc_kinds )
        bc_xmax = bc_kinds(j)
        if (bc_xmax == sll_p_periodic .and. (.not. pinfo % periodic)) cycle
        if (bc_xmax == sll_p_periodic .and. bc_xmin /= bc_xmax) cycle

        ! FIXME: spline 1D should be able to handle different BCs at xmin/xmax
        ! For now skip cases with bc_xmin /= bc_xmax
        if (bc_xmin /= bc_xmax) then
!          write(*,'(1i6)', advance='no') degree
!          write(*,'(2a6)', advance='no') bc_to_char( bc_xmin ), bc_to_char( bc_xmax )
!          write(*,'(a10,a8)') "---", "SKIP.."
          cycle
        end if

        ! Initialize test facility
        call test_facility % init( &
          profile_1d = profile_1d_cos, &
          degree  = degree , &
          ncells  = ncells , &
          bc_xmin = bc_xmin, &
!          bc_xmax = bc_xmax )
          bc_xmax = bc_xmax, &
          break_pts = [(pinfo%xmin+real(k,wp)*dx, k=0, ncells)] )

        ! Run tests
        call test_facility % evaluate_at_interpolation_points( max_norm_error )

        ! Calculate relative error norm from absolute one, check tolerance
        max_norm_error = max_norm_error / max_norm_profile
        success = (max_norm_error <= tol)

        ! Print test report to terminal on a single line
        write(*,'(1i6)', advance='no') degree
        write(*,'(2a6)', advance='no') bc_to_char( bc_xmin ), bc_to_char( bc_xmax )
        write(*,'(es10.1,a8)') max_norm_error, trim( success_to_string( success ))

        ! Free memory
        call test_facility % free()

        ! Update 'PASSED/FAILED' condition
        passed(1) = (passed(1) .and. success)

      end do  ! bc_xmax
    end do  ! bc_xmin
  end do  ! degree

  ! Deallocate local arrays
  deallocate( bc_kinds )

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 2: Spline should represent polynomial profiles exactly
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  write(*,*)
  write(*,'(a)') '---------------------------------------------------------------------------'
  write(*,'(a)') ' TEST 2: spline of order m exactly represents polynomial of same order     '
  write(*,'(a)') '---------------------------------------------------------------------------'
  write(*,*)

  ! Test tolerance: very small!
  tol      = 1e-14_wp
  tol_diff = 2e-13_wp ! larger for derivative

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(2) )
  bc_kinds(:) = [sll_p_hermite, sll_p_greville]

  ! Choose number of cells in grid
  ncells = 23

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = 1009

  ! Print constant parameters
  write(*,'(a,es10.1)') 'Relative error tolerance for f    : tol      = ', tol
  write(*,'(a,es10.1)') 'Relative error tolerance for ∂f/∂x: tol_diff = ', tol_diff
  write(*,'(a,i4)')     'Number of cells in grid: ncells = ', ncells
  write(*,'(a,i4)')     'Number of evaluation points (uniform grid): ', grid_dim
  write(*,*)
  write(*,'(a)') "Input:"
  write(*,'(a)') '  . deg   = spline degree'
  write(*,'(a)') '  . bcmin = boundary conditions at x=xmin [H|G]'
  write(*,'(a)') '  . bcmax = boundary conditions at x=xmax [H|G]'
  write(*,*)
  write(*,'(a)') "Output:"
  write(*,'(a)') '  .  err    = relative max-norm of error on f'
  write(*,'(a)') '  .  err_dx = relative max-norm of error on ∂f/∂x'
  write(*,'(a)') '  .  passed = "OK" if all errors <= tol, "FAIL" otherwise'
  write(*,*)
  write(*,'(a)') "Boundary conditions:"
  write(*,'(a)') '  .  H = Hermite'
  write(*,'(a)') '  .  G = Greville'
  write(*,*)

  ! Print table header
  write(*, '(3a6,2a10,a10)') "deg", "bcmin", "bcmax", "err", "err_dx", "passed"

  ! Create uniform grid of evaluation points
  call profile_1d_poly % init( 0 )  ! dummy, only need domain size
  call profile_1d_poly % get_info( pinfo )
  allocate( grid (grid_dim) )
  grid_dx = (pinfo%xmax-pinfo%xmin) / real( grid_dim-1, wp )
  grid    = [(pinfo%xmin + real(k,wp)*grid_dx, k=0, grid_dim-1)]

  ! Compute cell size in uniform grid
  dx = (pinfo%xmax-pinfo%xmin) / ncells

  ! Initialize 'PASSED/FAILED' condition
  passed(2) = .true.

  do degree = 1, 9

    ! Initialize polynomial profile with degree (deg1,deg2)
    call profile_1d_poly % init( degree )

    ! Estimate max-norm of profile (needed to compute relative error)
    max_norm_profile      = profile_1d_poly % max_norm()
    max_norm_profile_diff = profile_1d_poly % max_norm( diff=1 )

    do i = 1, size( bc_kinds )
      bc_xmin = bc_kinds(i)

      do j = 1, size( bc_kinds )
        bc_xmax = bc_kinds(j)

        ! FIXME: spline 1D should be able to handle different BCs at xmin/xmax
        ! For now skip cases with bc_xmin /= bc_xmax
        if (bc_xmin /= bc_xmax) then
!          write(*,'(1i6)', advance='no') degree
!          write(*,'(2a6)', advance='no') bc_to_char( bc_xmin ), bc_to_char( bc_xmax )
!          write(*,'(a10,a8)') "---", "SKIP.."
          cycle
        end if

        ! Initialize test facility
        call test_facility % init( &
          profile_1d = profile_1d_poly, &
          degree  = degree , &
          ncells  = ncells , &
          bc_xmin = bc_xmin, &
!          bc_xmax = bc_xmax )
          bc_xmax = bc_xmax, &
          break_pts = [(pinfo%xmin+real(k,wp)*dx, k=0, ncells)] )

        ! Run tests
        call test_facility % evaluate_on_1d_grid      ( grid, max_norm_error )
        call test_facility % evaluate_deriv_on_1d_grid( grid, max_norm_error_diff )

        ! Calculate relative error norms from absolute ones
        max_norm_error      = max_norm_error      / max_norm_profile
        max_norm_error_diff = max_norm_error_diff / max_norm_profile_diff

        ! Check tolerances
        success      = (max_norm_error      <= tol     )
        success_diff = (max_norm_error_diff <= tol_diff)

        ! Keep single success condition
        success = success .and. success_diff

        ! Print test report to terminal on a single line
        write(*,'(1i6)'    , advance='no') degree
        write(*,'(2a6)'    , advance='no') bc_to_char( bc_xmin ), bc_to_char( bc_xmax )
        write(*,'(2es10.1)', advance='no') max_norm_error, max_norm_error_diff
        write(*,'(a8)') trim( success_to_string( success ))

        ! Free memory
        call test_facility % free()

        ! Update 'PASSED/FAILED' condition
        passed(2) = (passed(2) .and. success)

      end do  ! bc_xmax
    end do  ! bc_xmin
  end do  ! degree

  ! Deallocate local arrays
  deallocate( bc_kinds )
  deallocate( grid     )

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 3: convergence analysis on cos*cos profile (with absolute error bound)
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  write(*,*)
  write(*,'(a)') '-----------------------------------------------------------------------------'
  write(*,'(a)') ' TEST 3: convergence analysis on cos profile (with absolute error bound)     '
  write(*,'(a)') '-----------------------------------------------------------------------------'
  write(*,*)

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(3) )
  bc_kinds(:) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  allocate( nx_list(5) )
  nx_list(:) = [10, 20, 40, 80, 160]

  ! Choose parameters of cos profile to be interpolated
  cos_n = 3
  cos_c = 0.3_wp

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = 20

  ! Print constant parameters
  write(*,'(a)') 'Profile: f(x) = cos( 2*pi * (n*x+c) )'
  write(*,'(a,i2)'    ) '  . n = ', cos_n
  write(*,'(a,f17.15)') '  . c = ', cos_c
  write(*,*)
  write(*,'(a,i4)')'Number of evaluation points (uniform grid): ', grid_dim
  write(*,*)
  write(*,'(a)') "Input:"
  write(*,'(a)') '  . deg   = spline degree (same along x1 and x2 directions)'
  write(*,'(a)') '  . bcmin = boundary conditions in x1 direction [P|H|G]'
  write(*,'(a)') '  . bcmax = boundary conditions in x2 direction [P|H|G]'
  write(*,'(a)') '  . nc    = number of grid cells'
  write(*,*)
  write(*,'(a)') "Output:"
  write(*,'(a)') '  . tol    = tolerance for f'
  write(*,'(a)') '  . err    = max-norm of error on f'
  write(*,'(a)') '  . tol_dx = tolerance for ∂f/∂x'
  write(*,'(a)') '  . err_dx = max-norm of error on ∂f/∂x'
  write(*,'(a)') '  . passed = "OK" if all errors <= tol, "FAIL" otherwise'
  write(*,*)
  write(*,'(a)') "Timing [s]:"
  write(*,'(a)') '  . t_init  = initialization of spline object'
  write(*,'(a)') '  . t_comp  = calculation of spline coefficients for interpolation'
  write(*,'(a)') '  . t_eval  = evaluation of function value S(x1,x2)'
  write(*,'(a)') '  . t_dx    = evaluation of x-derivative  ∂S(x1,x2)/∂x'
  write(*,*)
  write(*,'(a)') "Boundary conditions:"
  write(*,'(a)') '  . P = periodic'
  write(*,'(a)') '  . H = Hermite'
  write(*,'(a)') '  . G = Greville'
  write(*,*)

  ! Print table header
  write(*, '(a3,3a6,*(a10))') &
    "deg", "bcmin", "bcmax", "nc", &
    "tol", "err", "tol_dx1", "err_dx1", "passed", &
    "t_init", "t_comp", "t_eval", "t_dx"

  ! Initialize profile
  call profile_1d_cos % init( cos_n, cos_c )

  ! Extract information about 1D analytical profile
  call profile_1d_cos % get_info( pinfo )

  ! Estimate max-norm of profile (needed to compute relative error)
  max_norm_profile      = profile_1d_cos % max_norm()
  max_norm_profile_diff = profile_1d_cos % max_norm( diff=1 )

  ! Create uniform grid of evaluation points
  allocate( grid (grid_dim) )
  grid_dx = (pinfo%xmax-pinfo%xmin) / real( grid_dim-1, wp )
  do i = 1, grid_dim
    grid(i) = pinfo%xmin + real(i-1,wp)*grid_dx
  end do

  ! Initialize 'PASSED/FAILED' condition
  passed(3) = .true.

  do degree = 1, 9

    do i = 1, size( bc_kinds )
      bc_xmin = bc_kinds(i)

      do j = 1, size( bc_kinds )
        bc_xmax = bc_kinds(j)

        ! FIXME: spline 1D should be able to handle different BCs at xmin/xmax
        ! For now skip cases with bc_xmin /= bc_xmax
        if (bc_xmin /= bc_xmax) then
!          write(*,'(1i6)', advance='no') degree
!          write(*,'(2a6)', advance='no') bc_to_char( bc_xmin ), bc_to_char( bc_xmax )
!          write(*,'(a10,a8)') "---", "SKIP.."
          cycle
        end if

        do k = 1, size( nx_list )
          ncells = nx_list(k)

          ! TODO: remove this
          ! Compute cell size in uniform grid
          dx = (pinfo%xmax-pinfo%xmin) / ncells

          ! Initialize test facility
          call test_facility % init( &
            profile_1d = profile_1d_cos, &
            degree     = degree , &
            ncells     = ncells , &
            bc_xmin    = bc_xmin, &
!            bc_xmax    = bc_xmax )
            bc_xmax    = bc_xmax, &
            break_pts  = [(pinfo%xmin+real(k,wp)*dx, k=0, ncells)] )

          ! Determine error tolerances
          dx       = (pinfo%xmax - pinfo%xmin) / ncells
          tol      = sll_f_spline_1d_error_bound         ( profile_1d_cos, dx, degree )
          tol_diff = sll_f_spline_1d_error_bound_on_deriv( profile_1d_cos, dx, degree )

          ! Increase absolute tolerances if below machine precision
          ! NOTE: derivatives more sensitive to roundoff
          tol      = max( 1e-14_wp*max_norm_profile     , tol      )
          tol_diff = max( 1e-12_wp*max_norm_profile_diff, tol_diff )

          ! Run tests
          ! TODO: print numerical order of accuracy
          call test_facility % evaluate_on_1d_grid      ( grid, max_norm_error )
          call test_facility % evaluate_deriv_on_1d_grid( grid, max_norm_error_diff )

          ! Check tolerances
          success      = (max_norm_error      <= tol      )
          success_diff = (max_norm_error_diff <= tol_diff )

          ! Keep single success condition
          success = success .and. success_diff

          ! Print test report to terminal on a single line
          write(*,'(i3)'     , advance='no') degree
          write(*,'(2a6)'    , advance='no') bc_to_char( bc_xmin ), bc_to_char( bc_xmax )
          write(*,'(1i6)'    , advance='no') ncells
          write(*,'(2es10.1)', advance='no') tol     , max_norm_error
          write(*,'(2es10.1)', advance='no') tol_diff, max_norm_error_diff
          write(*,'(a10)'    , advance='no') success_to_string( success )
          write(*,'(4es10.1)') &
            test_facility % time_init               , &
            test_facility % time_compute_interpolant, &
            test_facility % time_eval_array         , &
            test_facility % time_eval_array_deriv

          ! Free memory
          call test_facility % free()

          ! Update 'PASSED/FAILED' condition
          passed(3) = (passed(3) .and. success)

        end do  ! nx1
        write(*,*)

      end do   ! bc2
    end do   ! bc1
  end do   ! deg1=deg2

  ! Deallocate local arrays
  deallocate( bc_kinds )
  deallocate( nx_list  )
  deallocate( grid     )

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  write(*,'(a)') '----------------'
!  write(*,'(a)') 'Test 0: '// success_to_string( passed(0) )
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

end program test_spline_1d
