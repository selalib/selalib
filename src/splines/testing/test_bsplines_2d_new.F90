program test_bsplines_2d_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: &
    f64

  use sll_m_constants, only: &
    sll_p_twopi

  use m_analytical_profiles_2d, only: &
    t_profile_2d_info, &
    c_analytical_profile_2d, &
    t_analytical_profile_2d_cos_cos, &
    t_analytical_profile_2d_poly

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite, &
    sll_p_greville

  use m_test_bsplines_2d, only: &
    t_bspline_2d_test_facility

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  ! Local variables
  type(t_profile_2d_info)               :: pinfo
  type(t_bspline_2d_test_facility)      :: test_facility
  type(t_analytical_profile_2d_cos_cos) :: profile_2d_cos_cos
  type(t_analytical_profile_2d_poly)    :: profile_2d_poly
  integer , allocatable                 :: bc_kinds(:)
  integer , allocatable                 :: nx_list(:)
  real(wp), allocatable                 :: grid_x1(:,:)
  real(wp), allocatable                 :: grid_x2(:,:)

  integer  :: deg1
  integer  :: deg2
  integer  :: nx1
  integer  :: nx2
  integer  :: bc1
  integer  :: bc2
  integer  :: i1, i2
  integer  :: grid_dim(2)
  real(wp) :: tol
  real(wp) :: tol_diff_x1
  real(wp) :: tol_diff_x2
  logical  :: equiv(3)
  logical  :: passed
  logical  :: success
  logical  :: success_diff_x1
  logical  :: success_diff_x2
  real(wp) :: max_norm_error        , max_norm_profile
  real(wp) :: max_norm_error_diff_x1, max_norm_profile_diff_x1
  real(wp) :: max_norm_error_diff_x2, max_norm_profile_diff_x2
  real(wp) :: dx1
  real(wp) :: dx2
  integer  :: j
  integer  :: cos_n1, cos_n2
  real(wp) :: cos_c1, cos_c2

  ! Initialize 'PASSED/FAILED' condition
  ! TODO: separate 'passed' value for each test and print to terminal
  passed = .true.

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 0: Check equivalence between scalar and array methods, and time them
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  write(*,*)
  write(*,'(a)') '---------------------------------------------------------------------------'
  write(*,'(a)') ' TEST 0: Check equivalence between scalar and array methods, and time them '
  write(*,'(a)') '---------------------------------------------------------------------------'
  write(*,*)

  ! Initialize profile: polynomial with some order
  call profile_2d_poly % init( 13, 12 )

  ! Choose some parameters
  nx1  = 19
  nx2  = 23
  deg1 = 5
  deg2 = 8
  bc1  = sll_p_greville
  bc2  = sll_p_hermite

  ! Initialize test facility
  call test_facility % init( &
    profile_2d = profile_2d_poly, &
    nx1        = nx1 , &
    nx2        = nx2 , &
    deg1       = deg1, &
    deg2       = deg2, &
    bc1        = bc1 , &
    bc2        = bc2  )

  ! Check equivalence between various methods
  call test_facility % check_equivalence_scalar_array_methods( equiv )

  write(*,'(a6,2a10,a10/)') 'method', 't_scalar', 't_array', 'equiv'

  write(*,'(a6,2es10.2,l8)') 'eval'      , &
    test_facility % time_eval            , &
    test_facility % time_eval_array      , equiv(1)

  write(*,'(a6,2es10.2,l8)') 'diff_x1'   , &
    test_facility % time_eval_diff1      , &
    test_facility % time_eval_diff1_array, equiv(2)

  write(*,'(a6,2es10.2,l8)') 'diff_x2'   , &
    test_facility % time_eval_diff2      , &
    test_facility % time_eval_diff2_array, equiv(3)

  ! Update 'PASSED/FAILED' condition
  passed = (passed .and. all( equiv ))

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 1: Evaluate spline at interpolation points (error should be zero)
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  write(*,*)
  write(*,'(a)') '-----------------------------------------------------------------------'
  write(*,'(a)') ' TEST 1: evaluate spline at interpolation points (error should be zero)'
  write(*,'(a)') '-----------------------------------------------------------------------'
  write(*,*)

  ! Test tolerance: very small!
  tol = 1e-14_wp

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(3) )
  bc_kinds(:) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  ! Choose number of knots in tensor grid
  nx1 = 10
  nx2 = 37

  ! Print constant parameters
  write(*,'(a)') 'Profile: f(x1,x2) = cos( 2*pi*x1 ) * cos( 2*pi*x2 )'
  write(*,*)
  write(*,'(a,e10.2)')  'Relative error tolerance: tol = ', tol
  write(*,'(a,i4)')     'Number of knots in grid : nx1 = ', nx1
  write(*,'(a,i4)')     '                          nx2 = ', nx2
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
  write(*, '(4a10,a12,a10)') "deg1", "deg2", "bc1", "bc2", "error", "passed"

  ! Initialize profile
  call profile_2d_cos_cos % init()

  ! Extract information about 2D analytical profile
  call profile_2d_cos_cos % get_info( pinfo )

  ! Estimate max-norm of profile (needed to compute relative error)
  max_norm_profile = profile_2d_cos_cos % max_norm()

  ! Cycle over spline degree
  do deg1 = 2, 9
    do deg2 = 2, 9

      ! Cycle over all kinds of boundary conditions in x1
      do i1 = 1, size( bc_kinds )
        bc1 = bc_kinds(i1)
        if (bc1 == sll_p_periodic .and. (.not. pinfo % x1_periodic)) cycle

        ! Cycle over all kinds of boundary conditions in x2
        do i2 = 1, size( bc_kinds )
          bc2 = bc_kinds(i2)
          if (bc2 == sll_p_periodic .and. (.not. pinfo % x2_periodic)) cycle

          ! Initialize test facility
          call test_facility % init( &
            profile_2d = profile_2d_cos_cos, &
            nx1        = nx1 , &
            nx2        = nx2 , &
            deg1       = deg1, &
            deg2       = deg2, &
            bc1        = bc1 , &
            bc2        = bc2  )

          ! Run tests
          call test_facility % evaluate_at_interpolation_points( max_norm_error )

          ! Calculate relative error norm from absolute one, check tolerance
          max_norm_error = max_norm_error / max_norm_profile
          success = (max_norm_error <= tol)

          ! Print test report to terminal on a single line
          write(*,'(2i10)', advance='no') deg1, deg2
          write(*,'(2a10)', advance='no') bc_to_char( bc1 ), bc_to_char( bc2 )
          write(*,'(e12.2,a8)') max_norm_error, trim( success_to_string( success ))

          ! Free memory
          call test_facility % free()

          ! Update 'PASSED/FAILED' condition
          passed = (passed .and. success)

        end do  ! bc2
      end do  ! bc1
    end do  ! deg2
  end do  ! deg1

  ! Deallocate local arrays
  deallocate( bc_kinds )

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

  ! Choose number of knots in tensor grid
  nx1 = 23
  nx2 = 18

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = [20, 20]

  ! Print constant parameters
  write(*,'(a,e10.2)')  'Relative error tolerance for f     : tol         = ', tol
  write(*,'(a,e10.2)')  'Relative error tolerance for ∂f/∂x1: tol_diff_x1 = ', tol_diff_x1
  write(*,'(a,e10.2)')  'Relative error tolerance for ∂f/∂x2: tol_diff_x2 = ', tol_diff_x2
  write(*,'(a,i4)')     'Number of knots in grid: nx1 = ', nx1
  write(*,'(a,i4)')     '                         nx2 = ', nx2
  write(*,'(a,i4,a,i4,a)')'Number of evaluation points (uniform grid): [', grid_dim(1), ',', grid_dim(2),']'
  write(*,*)
  write(*,'(a)') "Input:"
  write(*,'(a)') '  . deg1 = spline degree in x1 direction'
  write(*,'(a)') '  . deg2 = spline degree in x2 direction'
  write(*,'(a)') '  .  bc1 = boundary conditions in x1 direction [H|G]'
  write(*,'(a)') '  .  bc2 = boundary conditions in x2 direction [H|G]'
  write(*,*)
  write(*,'(a)') "Output:"
  write(*,'(a)') '  .  error     = relative max-norm of error on f'
  write(*,'(a)') '  .  error_dx1 = relative max-norm of error on ∂f/∂x1'
  write(*,'(a)') '  .  error_dx2 = relative max-norm of error on ∂f/∂x2'
  write(*,'(a)') '  .  passed    = "OK" if all errors <= tol, "FAIL" otherwise'
  write(*,*)
  write(*,'(a)') "Boundary conditions:"
  write(*,'(a)') '  .  H = Hermite'
  write(*,'(a)') '  .  G = Greville'
  write(*,*)

  ! Print table header
  write(*, '(4a10,3a12,a10)') "deg1", "deg2", "bc1", "bc2", &
    "error", "error_dx1", "error_dx2", "passed"

  ! Create uniform grid of evaluation points
  call profile_2d_poly % init( 0, 0 )  ! dummy, only need domain size
  call profile_2d_poly % get_info( pinfo )
  allocate( grid_x1 (grid_dim(1),grid_dim(2)) )
  allocate( grid_x2 (grid_dim(1),grid_dim(2)) )
  dx1 = (pinfo%x1_max-pinfo%x1_min) / real( grid_dim(1)-1, wp )
  dx2 = (pinfo%x2_max-pinfo%x2_min) / real( grid_dim(2)-1, wp )
  do i2 = 1, grid_dim(2)
    do i1 = 1, grid_dim(1)
      grid_x1(i1,i2) = pinfo%x1_min + real(i1-1,wp)*dx1
      grid_x2(i1,i2) = pinfo%x2_min + real(i2-1,wp)*dx2
    end do
  end do

  do deg1 = 2, 9
    do deg2 = 2, 9

      ! Initialize polynomial profile with degree (deg1,deg2)
      call profile_2d_poly % init( deg1, deg2 )

      ! Estimate max-norm of profile (needed to compute relative error)
      max_norm_profile         = profile_2d_poly % max_norm()
      max_norm_profile_diff_x1 = profile_2d_poly % max_norm( diff_x1=1 )
      max_norm_profile_diff_x2 = profile_2d_poly % max_norm( diff_x2=1 )

      do i1 = 1, size( bc_kinds )
        bc1 = bc_kinds(i1)

        do i2 = 1, size( bc_kinds )
          bc2 = bc_kinds(i2)

          ! Initialize test facility
          call test_facility % init( &
            profile_2d = profile_2d_poly, &
            nx1        = nx1 , &
            nx2        = nx2 , &
            deg1       = deg1, &
            deg2       = deg2, &
            bc1        = bc1 , &
            bc2        = bc2 )

          ! Run tests
          call test_facility % evaluate_on_2d_grid( grid_x1, grid_x2, max_norm_error )
          call test_facility % evaluate_grad_on_2d_grid( grid_x1, grid_x2, &
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
          write(*,'(6i10)'  , advance='no') deg1, deg2
          write(*,'(2a10)'  , advance='no') bc_to_char( bc1 ), bc_to_char( bc2 )
          write(*,'(3e12.2)', advance='no') &
            max_norm_error, max_norm_error_diff_x1, max_norm_error_diff_x2
          write(*,'(a8)') trim( success_to_string( success ))

          ! Free memory
          call test_facility % free()

          ! Update 'PASSED/FAILED' condition
          passed = (passed .and. success)

        end do  ! bc2
      end do  ! bc1
    end do  ! deg2
  end do  ! deg1

  ! Deallocate local arrays
  deallocate( bc_kinds )
  deallocate( grid_x1  )
  deallocate( grid_x2  )

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
  nx_list(:) = [11, 21, 41, 81, 161]

  ! Choose parameters of cos*cos profile to be interpolated
  cos_n1 = 3
  cos_n2 = 3
  cos_c1 = 0.3_wp
  cos_c2 = 0.7_wp

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = [20,20]

  ! Print constant parameters
  write(*,'(a)') 'Profile: f(x1,x2) = cos( 2*pi * (n1*x1+c1) ) * cos( 2*pi * (n2*x2+c2) )'
  write(*,'(a,i2)') '  . n1 = ', cos_n1
  write(*,'(a,i2)') '  . n2 = ', cos_n2
  write(*,'(a,f17.15)')  '  . c1 = ', cos_c1
  write(*,'(a,f17.15)')  '  . c2 = ', cos_c2
  write(*,*)
  write(*,'(a,i4,a,i4,a)')'Number of evaluation points (uniform grid): [', grid_dim(1), ',', grid_dim(2),']'
  write(*,*)
  write(*,'(a)') "Input:"
  write(*,'(a)') '  . deg = spline degree (same along x1 and x2 directions)'
  write(*,'(a)') '  . bc1 = boundary conditions in x1 direction [P|H|G]'
  write(*,'(a)') '  . bc2 = boundary conditions in x2 direction [P|H|G]'
  write(*,'(a)') '  . nx  = number of grid points along each direction (nx1=nx2)'
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

  ! Create uniform grid of evaluation points
  allocate( grid_x1 (grid_dim(1),grid_dim(2)) )
  allocate( grid_x2 (grid_dim(1),grid_dim(2)) )
  dx1 = (pinfo%x1_max-pinfo%x1_min) / real( grid_dim(1)-1, wp )
  dx2 = (pinfo%x2_max-pinfo%x2_min) / real( grid_dim(2)-1, wp )
  do i2 = 1, grid_dim(2)
    do i1 = 1, grid_dim(1)
      grid_x1(i1,i2) = pinfo%x1_min + real(i1-1,wp)*dx1
      grid_x2(i1,i2) = pinfo%x2_min + real(i2-1,wp)*dx2
    end do
  end do

  do deg1 = 2, 9

    ! Use same spline degree along both directions
    deg2 = deg1

    do i1 = 1, size( bc_kinds )
      bc1 = bc_kinds(i1)

      do i2 = 1, size( bc_kinds )
        bc2 = bc_kinds(i2)

        do j = 1, size( nx_list )
          nx1 = nx_list(j)

          ! Use same number of knots along both directions
          nx2 = nx1

          ! Initialize test facility
          call test_facility % init( &
            profile_2d = profile_2d_cos_cos, &
            nx1        = nx1 , &
            nx2        = nx2 , &
            deg1       = deg1, &
            deg2       = deg2, &
            bc1        = bc1 , &
            bc2        = bc2 )

          ! Determine error tolerances
          dx1 = (pinfo%x1_max - pinfo%x1_min) / nx1
          dx2 = (pinfo%x2_max - pinfo%x2_min) / nx2
          tol         = error_bound        ( profile_2d_cos_cos, dx1, dx2, deg1, deg2 )
          tol_diff_x1 = error_bound_diff_x1( profile_2d_cos_cos, dx1, dx2, deg1, deg2 )
          tol_diff_x2 = error_bound_diff_x2( profile_2d_cos_cos, dx1, dx2, deg1, deg2 )

          ! Increase absolute tolerances if below machine precision
          ! NOTE: derivatives more sensitive to roundoff
          tol         = max( 1e-14_wp*max_norm_profile        , tol         )
          tol_diff_x1 = max( 1e-12_wp*max_norm_profile_diff_x1, tol_diff_x1 )
          tol_diff_x2 = max( 1e-12_wp*max_norm_profile_diff_x2, tol_diff_x2 )

          ! Run tests
          ! TODO: print numerical order of accuracy
          call test_facility % evaluate_on_2d_grid( grid_x1, grid_x2, max_norm_error )
          call test_facility % evaluate_grad_on_2d_grid( grid_x1, grid_x2, &
            max_norm_error_diff_x1, max_norm_error_diff_x2 )

          ! Check tolerances
          success         = (max_norm_error         <= tol        )
          success_diff_x1 = (max_norm_error_diff_x1 <= tol_diff_x1)
          success_diff_x2 = (max_norm_error_diff_x2 <= tol_diff_x2)

          ! Keep single success condition
          success = success .and. success_diff_x1 .and. success_diff_x2

          ! Print test report to terminal on a single line
          write(*,'(i3)'     , advance='no') deg1
          write(*,'(2a6)'    , advance='no') bc_to_char( bc1 ), bc_to_char( bc2 )
          write(*,'(2i6)'    , advance='no') nx1
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
          passed = (passed .and. success)

        end do  ! nx1
        write(*,*)

      end do   ! bc2
    end do   ! bc1
  end do   ! deg1=deg2

  ! Deallocate local arrays
  deallocate( bc_kinds )
  deallocate( nx_list  )
  deallocate( grid_x1  )
  deallocate( grid_x2  )

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  ! CTest key
  if (passed) then
    write(*,*) "PASSED"
  else
    write(*,*) "FAILED"
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
  ! Error bound in max norm for spline interpolation of periodic functions from:
  !
  ! V M Tihomirov 1969 Math. USSR Sb. 9 275
  ! https://doi.org/10.1070/SM1969v009n02ABEH002052 (page 286, bottom)
  !
  ! Yu. S. Volkov and Yu. N. Subbotin
  ! https://doi.org/10.1134/S0081543815020236 (equation 14)
  !
  ! Also applicable to first derivative by passing deg-1 instead of deg
  ! Volkov & Subbotin 2015, eq. 15
  !-----------------------------------------------------------------------------
  pure function sll_f_spline_1d_error_bound( h, deg, norm_f ) result( norm_e )
    real(wp), intent(in) :: h
    integer , intent(in) :: deg
    real(wp), intent(in) :: norm_f
    real(wp) :: norm_e

    real(wp), parameter :: k(1:10) = &
         [     1.0_wp /          2.0_wp, &
               1.0_wp /          8.0_wp, &
               1.0_wp /         24.0_wp, &
               5.0_wp /        384.0_wp, &
               1.0_wp /        240.0_wp, &
              61.0_wp /      46080.0_wp, &
              17.0_wp /      40320.0_wp, &
             277.0_wp /    2064384.0_wp, &
              31.0_wp /     725760.0_wp, &
           50521.0_wp / 3715891200.0_wp ]

    norm_e = k(deg+1) * h**(deg+1) * norm_f

  end function sll_f_spline_1d_error_bound

  !-----------------------------------------------------------------------------
  function error_bound( profile_2d, dx1, dx2, deg1, deg2 ) result( max_error )
    class(c_analytical_profile_2d), intent(in) :: profile_2d
    real(wp)                      , intent(in) :: dx1
    real(wp)                      , intent(in) :: dx2
    integer                       , intent(in) :: deg1
    integer                       , intent(in) :: deg2
    real(wp) :: max_error

    real(wp) :: max_norm1
    real(wp) :: max_norm2

    max_norm1 = profile_2d % max_norm( deg1+1, 0      )
    max_norm2 = profile_2d % max_norm( 0     , deg2+1 )

    max_error = sll_f_spline_1d_error_bound( dx1, deg1, max_norm1 ) &
              + sll_f_spline_1d_error_bound( dx2, deg2, max_norm2 )

  end function error_bound

  !-----------------------------------------------------------------------------
  function error_bound_diff_x1( profile_2d, dx1, dx2, deg1, deg2 ) result( max_error )
    class(c_analytical_profile_2d), intent(in) :: profile_2d
    real(wp)                      , intent(in) :: dx1
    real(wp)                      , intent(in) :: dx2
    integer                       , intent(in) :: deg1
    integer                       , intent(in) :: deg2
    real(wp) :: max_error

    real(wp) :: max_norm1
    real(wp) :: max_norm2

    max_norm1 = profile_2d % max_norm( deg1+1, 0      )
    max_norm2 = profile_2d % max_norm( 0     , deg2+1 )

    max_error = sll_f_spline_1d_error_bound( dx1, deg1-1, max_norm1 ) &
              + sll_f_spline_1d_error_bound( dx2, deg2  , max_norm2 )

  end function error_bound_diff_x1

  !-----------------------------------------------------------------------------
  function error_bound_diff_x2( profile_2d, dx1, dx2, deg1, deg2 ) result( max_error )
    class(c_analytical_profile_2d), intent(in) :: profile_2d
    real(wp)                      , intent(in) :: dx1
    real(wp)                      , intent(in) :: dx2
    integer                       , intent(in) :: deg1
    integer                       , intent(in) :: deg2
    real(wp) :: max_error

    real(wp) :: max_norm1
    real(wp) :: max_norm2

    max_norm1 = profile_2d % max_norm( deg1+1, 0      )
    max_norm2 = profile_2d % max_norm( 0     , deg2+1 )

    max_error = sll_f_spline_1d_error_bound( dx1, deg1  , max_norm1 ) &
              + sll_f_spline_1d_error_bound( dx2, deg2-1, max_norm2 )

  end function error_bound_diff_x2

end program test_bsplines_2d_new
