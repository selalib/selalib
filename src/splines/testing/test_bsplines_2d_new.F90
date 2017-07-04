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
  logical  :: passed
  logical  :: success
  real(wp) :: max_norm_error
  real(wp) :: max_norm_profile
  real(wp) :: dx1
  real(wp) :: dx2
  integer  :: j

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 1: Evaluate spline at interpolation points (error should be zero)
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  ! TODO: print test explanation to terminal

  ! Test tolerance: very small!
  tol = 1e-14_wp

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(3) )
  bc_kinds(:) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  ! Initialize 'PASSED/FAILED' condition
  passed = .true.

  ! Print report header
  call test_facility % print_header()

  ! Initialize profile
  call profile_2d_cos_cos % init()

  ! Extract information about 2D analytical profile
  call profile_2d_cos_cos % get_info( pinfo )

  ! Estimate max-norm of profile (needed to compute relative error)
  max_norm_profile = profile_2d_cos_cos % max_norm()

  ! Choose number of knots in tensor grid
  nx1 = 10
  nx2 = 37

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
          success = (max_norm_error <= tol * max_norm_profile)

          ! Print test report to terminal on a single line
          write(*,'(6i10)', advance='no') nx1, nx2, deg1, deg2, bc1, bc2
          write(*,'(e12.2,L12)') max_norm_error / max_norm_profile, success

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

  ! TODO: print test explanation to terminal
  write(*,*)
  write(*,*)

  ! Test tolerance: very small!
  tol = 1e-14_wp

  ! Print report header
  call test_facility % print_header()

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(2) )
  bc_kinds(:) = [sll_p_hermite, sll_p_greville]

  ! Choose number of knots in tensor grid
  nx1 = 23
  nx2 = 18

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = [20, 20]

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
      max_norm_profile = profile_2d_poly % max_norm()

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
          success = (max_norm_error <= tol * max_norm_profile)

          ! Print test report to terminal on a single line
          write(*,'(6i10)', advance='no') nx1, nx2, deg1, deg2, bc1, bc2
          write(*,'(e12.2,L12)') max_norm_error / max_norm_profile, success

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
  ! TEST 3: convergence analysis on cos*cos profile (with error bound)
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  ! TODO: print test explanation to terminal
  write(*,*)
  write(*,*)

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(3) )
  bc_kinds(:) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  allocate( nx_list(5) )
  nx_list(:) = [11, 21, 41, 81, 161]

  ! Print report header
  call test_facility % print_header()

  ! Initialize profile
  call profile_2d_cos_cos % init( n1=3, n2=3, c1=0.3_wp, c2=0.7_wp )

  ! Extract information about 2D analytical profile
  call profile_2d_cos_cos % get_info( pinfo )

  ! Estimate max-norm of profile (needed to compute relative error)
  max_norm_profile = profile_2d_cos_cos % max_norm()

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = [20,20]

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

          ! Determine error tolerance
          dx1 = (pinfo%x1_max - pinfo%x1_min) / nx1
          dx2 = (pinfo%x2_max - pinfo%x2_min) / nx2
          tol = error_bound( profile_2d_cos_cos, dx1, dx2, deg1, deg2 )
          tol = max( 1e-14_wp, tol )
          !write(*,*) tol

          ! Run tests
          ! TODO: print tolerance
          ! TODO: print numerical order of accuracy
          call test_facility % evaluate_on_2d_grid( grid_x1, grid_x2, max_norm_error )
          success = (max_norm_error <= tol * max_norm_profile)

          ! Print test report to terminal on a single line
          write(*,'(6i10)', advance='no') nx1, nx2, deg1, deg2, bc1, bc2
          write(*,'(e12.2,L12)') max_norm_error / max_norm_profile, success

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
    real(wp), intent(in   ) :: dx1
    real(wp), intent(in   ) :: dx2
    integer , intent(in   ) :: deg1
    integer , intent(in   ) :: deg2
    real(wp) :: max_error

    real(wp) :: max_norm1
    real(wp) :: max_norm2

    max_norm1 = profile_2d % max_norm( deg1+1, 0      )
    max_norm2 = profile_2d % max_norm( 0     , deg2+1 )

    max_error = sll_f_spline_1d_error_bound( dx1, deg1, max_norm1 ) &
              + sll_f_spline_1d_error_bound( dx2, deg2, max_norm2 )

  end function error_bound

end program test_bsplines_2d_new
