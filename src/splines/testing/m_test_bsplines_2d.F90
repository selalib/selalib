module m_test_bsplines_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: &
    f64

  use m_analytical_profiles_2d_base, only: &
    t_profile_2d_info, &
    c_analytical_profile_2d

  use sll_m_boundary_condition_descriptors, only: &
       sll_p_periodic, &
       sll_p_hermite, &
       sll_p_greville

  use sll_m_bspline_2d, only: &
     sll_t_bspline_2d,        &
     sll_s_bspline_2d_init,   &
     sll_s_bspline_2d_free,   &
     sll_s_compute_bspline_2d,                  &
     sll_f_interpolate_value_2d,                &
     sll_s_interpolate_array_values_2d,         &
     sll_f_interpolate_derivative_x1_2d,        &
     sll_f_interpolate_derivative_x2_2d,        &
     sll_s_interpolate_array_derivatives_x1_2d, &
     sll_s_interpolate_array_derivatives_x2_2d


  implicit none

  public :: &
    s_check_all_bcs

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> List of available boundary conditions
  integer, parameter :: bc_kinds(3) = &
                                 [sll_p_periodic, sll_p_hermite, sll_p_greville]

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_check_all_bcs( profile_2d, degree1, degree2, nx1, nx2 )
    class( c_analytical_profile_2d ), intent(in) :: profile_2d
    integer                         , intent(in) :: degree1
    integer                         , intent(in) :: degree2
    integer                         , intent(in) :: nx1
    integer                         , intent(in) :: nx2

    integer :: bc1
    integer :: bc2
    integer :: i1
    integer :: i2
    type(t_profile_2d_info) :: pinfo
    type(sll_t_bspline_2d ) :: bspline_2d

    ! Extract information about 2D analytical profile
    call profile_2d % get_info( pinfo )

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

        ! Initialize spline, each time with new boundary conditions in x1 and x2
        ! TODO: Hermite BCs require additional input data
        call sll_s_bspline_2d_init( bspline_2d, &
          nx1, &
          nx2, &
          degree1, &
          degree2, &
          pinfo % x1_min, &
          pinfo % x2_min, &
          pinfo % x1_max, &
          pinfo % x2_max, &
          bc1, &
          bc2 )

        write(*,*) "2D B-spline initialized with parameters:"
        write(*,*) "  . nx1     = ", nx1
        write(*,*) "  . nx2     = ", nx2
        write(*,*) "  . degree1 = ", degree1
        write(*,*) "  . degree2 = ", degree2
        write(*,*) "  . x1_min  = ", pinfo % x1_min
        write(*,*) "  . x2_min  = ", pinfo % x2_min
        write(*,*) "  . x1_max  = ", pinfo % x1_max
        write(*,*) "  . x2_max  = ", pinfo % x2_max
        write(*,*) "  . bc1     = ", bc1
        write(*,*) "  . bc2     = ", bc2
        write(*,*)

        ! Interpolate analytical 2D profile and run some tests
        call s_interpolate_and_evaluate( profile_2d, bspline_2d )

        ! Free spline memory
        call sll_s_bspline_2d_free( bspline_2d )

      end do
    end do

  end subroutine s_check_all_bcs

  !-----------------------------------------------------------------------------
  subroutine s_interpolate_and_evaluate( profile_2d, bspline_2d )
    class( c_analytical_profile_2d ), intent(in   ) :: profile_2d
    type ( sll_t_bspline_2d        ), intent(inout) :: bspline_2d

    type(t_profile_2d_info) :: pinfo
    real(wp), pointer       :: tau1(:)
    real(wp), pointer       :: tau2(:)
    real(wp), allocatable   :: gtau(:,:)  ! Profile values at interp. points
!    real(wp), allocatable   :: stau(:,:)  ! Spline  values at interp. points
    integer                 :: nipts1
    integer                 :: nipts2
    integer                 :: deg1
    integer                 :: deg2
    integer                 :: i1, j1, s1
    integer                 :: i2, j2, s2
    real(wp)                :: error
    real(wp)                :: max_norm_error

    real(wp), allocatable :: bc1_min(:,:)
    real(wp), allocatable :: bc1_max(:,:)
    real(wp), allocatable :: bc2_min(:,:)
    real(wp), allocatable :: bc2_max(:,:)

    ! Get profile info
    call profile_2d % get_info( pinfo )

    ! Get spline degree
    deg1 = bspline_2d % bs1 % deg
    deg2 = bspline_2d % bs2 % deg

    ! Get spline interpolation points
    tau1 => bspline_2d % bs1 % tau
    tau2 => bspline_2d % bs2 % tau

    ! Store number of interpolation points
    nipts1 = size( tau1 )
    nipts2 = size( tau2 )
    print *, "nipts1 = ", nipts1
    print *, "nipts2 = ", nipts2

    allocate( gtau (nipts1,nipts2) )
!    allocate( stau (nipts1,nipts2) )

    ! Evaluate analytical profile at interpolation points
    do i2 = 1, nipts2
      do i1 = 1, nipts1
        gtau(i1,i2) = profile_2d % eval( tau1(i1), tau2(i2) )
      end do
    end do

    ! If needed, evaluate derivatives at domain boundary
    ! TODO: add cross-derivatives at corners
    if (bspline_2d % bs1 % bc_type == sll_p_hermite) then
      allocate( bc1_min (deg1/2, nipts2) )
      allocate( bc1_max (deg1/2, nipts2) )
      s1 = 1-modulo(deg1,2) ! shift = 1 for even order, 0 for odd order
      do i2 = 1, nipts2
        do j1 = 1, deg1/2
          bc1_min(j1,i2) = profile_2d % eval( pinfo%x1_min, tau2(i2), diff_x1=j1-s1 )
          bc1_max(j1,i2) = profile_2d % eval( pinfo%x1_max, tau2(i2), diff_x1=j1-s1 )
        end do
      end do
    end if

    if (bspline_2d % bs2 % bc_type == sll_p_hermite) then
      allocate( bc2_min (deg2/2, nipts1) )
      allocate( bc2_max (deg2/2, nipts1) )
      s2 = 1-modulo(deg2,2) ! shift = 1 for even order, 0 for odd order
      do i1 = 1, nipts1
        do j2 = 1, deg2/2
          bc2_min(j2,i1) = profile_2d % eval( tau1(i1), pinfo%x2_min, diff_x2=j2-s2 )
          bc2_max(j2,i1) = profile_2d % eval( tau1(i1), pinfo%x2_max, diff_x2=j2-s2 )
        end do
      end do
    end if

    ! Compute 2D spline that interpolates analytical 2D profile at points above
    ! TODO: for Hermite we should also pass val1[2]_min[max]
    call sll_s_compute_bspline_2d( bspline_2d, gtau )
!    call sll_s_compute_bspline_2d( bspline_2d, gtau, &
!      val1_min = bc1_min, &
!      val1_max = bc1_max, &
!      val2_min = bc2_min, &
!      val2_max = bc2_max )

    ! Evaluate 2D spline at interpolation points: error should be zero
    ! TODO: max index correct only for cubic splines
    max_norm_error = 0.0_wp
    do i2 = 1, nipts2-2
      do i1 = 1, nipts1-2
        error = gtau(i1,i2) &
              - sll_f_interpolate_value_2d( bspline_2d, tau1(i1), tau2(i2) )
        max_norm_error = max( max_norm_error, abs( error ) )
      end do
    end do

    write(*,*) "Evaluate 2D spline at interpolation points: error should be zero"
    write(*,*) "max_norm_error = ", max_norm_error
    write(*,*)

    ! Deallocate arrays
    deallocate( gtau )
    if (allocated(bc1_min)) deallocate( bc1_min )
    if (allocated(bc1_max)) deallocate( bc1_max )
    if (allocated(bc2_min)) deallocate( bc2_min )
    if (allocated(bc2_max)) deallocate( bc2_max )

  end subroutine s_interpolate_and_evaluate

end module m_test_bsplines_2d
