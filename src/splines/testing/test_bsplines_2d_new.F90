program test_bsplines_2d_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: &
    f64

  use m_analytical_profiles_2d, only: &
    t_profile_2d_info, &
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
  integer, allocatable                  :: bc_kinds(:)

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

  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  ! TEST 1: Evaluate spline at interpolation points (error should be zero)
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  ! TODO: print test explanation to terminal

  ! Test tolerance: very small!
  tol = 1e-14_wp

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(3) )
  bc_kinds(1:3) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  ! Initialize 'PASSED/FAILED' condition
  passed = .true.

  ! Print report header
  call test_facility % print_header()

  ! Initialize profile
  call profile_2d_cos_cos % init()

  ! Extract information about 2D analytical profile
  call profile_2d_cos_cos % get_info( pinfo )

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
          call test_facility % evaluate_at_interpolation_points( tol, success )

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

  ! Test tolerance: small but not extremely so
  tol = 1e-13_wp

  ! Print report header
  call test_facility % print_header()

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(2) )
  bc_kinds(1:2) = [sll_p_hermite, sll_p_greville]

  ! Choose number of knots in tensor grid
  nx1 = 23
  nx2 = 18

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = [20, 20]

  do deg1 = 2, 9
    do deg2 = 2, 9

      ! Initialize polynomial profile with degree (deg1,deg2)
      call profile_2d_poly % init( deg1, deg2 )

      do i1 = 1, size( bc_kinds )
        bc1 = bc_kinds(i1)

        do i2 = 1, size( bc_kinds )
          bc2 = bc_kinds(i2)

          call test_facility % init( &
            profile_2d = profile_2d_poly, &
            nx1        = nx1 , &
            nx2        = nx2 , &
            deg1       = deg1, &
            deg2       = deg2, &
            bc1        = bc1 , &
            bc2        = bc2 )

          call test_facility % evaluate_func_and_grad_on_uniform_grid( &
            grid_dim, tol, success )

          call test_facility % free()

          passed = (passed .and. success)

        end do  ! bc2
      end do  ! bc1
    end do  ! deg2
  end do  ! deg1

  ! Deallocate local arrays
  deallocate( bc_kinds )
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  ! CTest key
  if (passed) then
    write(*,*) "PASSED"
  else
    write(*,*) "FAILED"
  end if

end program test_bsplines_2d_new
