program test_bsplines_2d_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: &
    f64

  use sll_m_constants, only: &
    sll_p_twopi

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
  integer, allocatable                  :: nx_list(:)

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
  tol = 2e-13_wp

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

  do deg1 = 2, 9
    do deg2 = 2, 9

      ! Initialize polynomial profile with degree (deg1,deg2)
      call profile_2d_poly % init( deg1, deg2 )

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
          call test_facility % evaluate_func_and_grad_on_uniform_grid( &
            grid_dim, tol, success )

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
  ! TEST 3: convergence analysis on cos*cos profile (with error bound)
  ! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  ! TODO: print test explanation to terminal
  write(*,*)
  write(*,*)

  ! Choose boundary conditions to be tested
  allocate( bc_kinds(3) )
  bc_kinds(:) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  allocate( nx_list(4) )
  nx_list(:) = [13, 26, 53, 78]

  ! Print report header
  call test_facility % print_header()

  ! Initialize profile
  call profile_2d_cos_cos % init()

  ! Extract information about 2D analytical profile
  call profile_2d_cos_cos % get_info( pinfo )

  ! Choose dimension of uniform grid of evaluation points
  grid_dim = [20,20]

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
          tol = error_bound_cos_cos( dx1, dx2, deg1, deg2 )
          tol = max( 1e-14_wp, tol )

          ! Run tests
          ! TODO: print tolerance
          ! TODO: evaluate all grid points at once
          ! TODO: print numerical order of accuracy
          call test_facility % evaluate_func_and_grad_on_uniform_grid( &
            grid_dim, tol, success )

          ! Free memory
          call test_facility % free()

          ! Update 'PASSED/FAILED' condition
          passed = (passed .and. success)

        end do
        write(*,*)

      end do
    end do
  end do

  ! Deallocate local arrays
  deallocate( bc_kinds )
  deallocate( nx_list  )

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
  pure function taylor_error_term_cos_cos( k1, k2 ) result( e )
    integer, intent(in) :: k1
    integer, intent(in) :: k2
    real(wp) :: e

    e = sll_p_twopi**(k1+k2) * dx1**k1 * dx2**k2 &
                / (factorial( k1 ) * factorial( k2 ))

  end function taylor_error_term_cos_cos

  !-----------------------------------------------------------------------------
  pure function error_bound_cos_cos( dx1, dx2, deg1, deg2 ) result( max_error )
    real(wp), intent(in   ) :: dx1
    real(wp), intent(in   ) :: dx2
    integer , intent(in   ) :: deg1
    integer , intent(in   ) :: deg2
    real(wp) :: max_error

!    integer :: k1, k2
!
!    max_error = 0.0_wp
!
!    k2 = deg2+1
!    do k1 = 0, deg1
!      max_error = max_error + taylor_error_term_cos_cos( k1, k2 )
!    end do
!
!    k1 = deg1+1
!    do k2 = 0, deg2
!      max_error = max_error + taylor_error_term_cos_cos( k1, k2 )
!    end do
!
!    max_error = max_error + taylor_error_term_cos_cos( deg1+1, deg2+1 )

    max_error = taylor_error_term_cos_cos( deg1+1, 0 ) &
              + taylor_error_term_cos_cos( 0, deg2+1 )

  end function error_bound_cos_cos

end program test_bsplines_2d_new
