program test_bsplines_2d_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: &
    f64

  use m_analytical_profiles_2d, only: &
    t_profile_2d_info, &
    t_analytical_profile_2d_cos_cos

  use sll_m_boundary_condition_descriptors, only: &
       sll_p_periodic, &
       sll_p_hermite, &
       sll_p_greville

  use m_test_bsplines_2d, only: &
    t_bspline_2d_test_facility

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! List of available boundary conditions
  integer, parameter :: bc_kinds(3) = &
                                 [sll_p_periodic, sll_p_hermite, sll_p_greville]

  ! Local variables
  type(t_profile_2d_info)               :: pinfo
  type(t_analytical_profile_2d_cos_cos) :: profile_2d_cos_cos
  type(t_bspline_2d_test_facility)      :: test_facility
  integer :: deg1
  integer :: deg2
  integer :: nx1
  integer :: nx2
  integer :: bc1
  integer :: bc2
  integer :: i1, i2
  logical :: passed
  logical :: successful

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
        if (bc1 == sll_p_periodic .and. (.not. pinfo % x1_periodic)) then
          cycle
        end if

        ! Cycle over all kinds of boundary conditions in x2
        do i2 = 1, size( bc_kinds )
          bc2 = bc_kinds(i2)
          if (bc2 == sll_p_periodic .and. (.not. pinfo % x2_periodic)) then
            cycle
          end if

          ! Initialize test facility
          call test_facility % init( &
            profile_2d = profile_2d_cos_cos, &
            nx1        = nx1, &
            nx2        = nx2, &
            deg1       = deg1, &
            deg2       = deg2, &
            bc1        = bc1, &
            bc2        = bc2 )

          ! Run tests
          call test_facility % run_tests( verbose=.false., successful=successful )

          passed = (passed .and. successful)

        end do  ! bc2
      end do  ! bc1
    end do  ! deg2
  end do  ! deg1

  ! CTest key
  if (passed) then
    write(*,*) "PASSED"
  else
    write(*,*) "FAILED"
  end if

end program test_bsplines_2d_new
