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
    t_bspline_2d_test_facility

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Test report
  type :: t_test_report
    character(len=256) :: name
    real(wp)           :: tol
    real(wp)           :: error
    logical            :: passed
  end type t_test_report


  !> Type for running test
  type :: t_bspline_2d_test_facility

    class(c_analytical_profile_2d), pointer :: profile_2d
    integer                                 :: nx1
    integer                                 :: nx2
    integer                                 :: deg1
    integer                                 :: deg2
    integer                                 :: bc1
    integer                                 :: bc2
    type(t_test_report)                     :: report

    type(sll_t_bspline_2d)  :: bspline_2d
    real(wp), allocatable   :: gtau(:,:)  ! Profile values at interp. points

  contains

    procedure, nopass :: print_header
    procedure         :: init
    procedure         :: run_tests

  end type t_bspline_2d_test_facility

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  subroutine print_header()

    integer, parameter :: n_tests = 1
    integer            :: i
    character(len=6)   :: is

    write(*, '(6a10)', advance='no') &
      "nx1", "nx2", "degree1", "degree2", "bc1", "bc2"

    do i = 1, n_tests
      write(is,'(i1)') i
      write( *,'(a12,a12)', advance='no') &
        't'//trim(is)//'%error', &
        't'//trim(is)//'%passed'
    end do
    write(*,*)

  end subroutine print_header

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
    real(wp), pointer       :: tau1(:)
    real(wp), pointer       :: tau2(:)
    integer                 :: nipts1
    integer                 :: nipts2
    integer                 :: i1, j1, s1
    integer                 :: i2, j2, s2

    real(wp), allocatable :: bc1_min(:,:)
    real(wp), allocatable :: bc1_max(:,:)
    real(wp), allocatable :: bc2_min(:,:)
    real(wp), allocatable :: bc2_max(:,:)
    real(wp), allocatable :: val_corners(:,:,:)

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
      self % bc1, &
      self % bc2 )

    ! Get spline interpolation points
    tau1 => self % bspline_2d % bs1 % tau
    tau2 => self % bspline_2d % bs2 % tau

    ! Store number of interpolation points
    nipts1 = size( tau1 )
    nipts2 = size( tau2 )
!    print *, "nipts1 = ", nipts1
!    print *, "nipts2 = ", nipts2

    ! Evaluate analytical profile at interpolation points
    allocate( self % gtau (nipts1,nipts2) )
    do i2 = 1, nipts2
      do i1 = 1, nipts1
        self % gtau(i1,i2) = self % profile_2d % eval( tau1(i1), tau2(i2) )
      end do
    end do

    ! If needed, evaluate x1 derivatives at (x1_min,x2) and (x1_max,x2)
    if (self % bspline_2d % bs1 % bc_type == sll_p_hermite) then
      allocate( bc1_min (deg1/2, nipts2) )
      allocate( bc1_max (deg1/2, nipts2) )
      s1 = 1-modulo(deg1,2) ! shift = 1 for even order, 0 for odd order
      do i2 = 1, nipts2
        do j1 = 1, deg1/2
          bc1_min(j1,i2) = self % profile_2d % eval( info%x1_min, tau2(i2), diff_x1=j1-s1 )
          bc1_max(j1,i2) = self % profile_2d % eval( info%x1_max, tau2(i2), diff_x1=j1-s1 )
        end do
      end do
    end if

    ! If needed, evaluate x2 derivatives at (x1,x2_min) and (x1,x2_max)
    if (self % bspline_2d % bs2 % bc_type == sll_p_hermite) then
      allocate( bc2_min (deg2/2, nipts1) )
      allocate( bc2_max (deg2/2, nipts1) )
      s2 = 1-modulo(deg2,2) ! shift = 1 for even order, 0 for odd order
      do i1 = 1, nipts1
        do j2 = 1, deg2/2
          bc2_min(j2,i1) = self % profile_2d % eval( tau1(i1), info%x2_min, diff_x2=j2-s2 )
          bc2_max(j2,i1) = self % profile_2d % eval( tau1(i1), info%x2_max, diff_x2=j2-s2 )
        end do
      end do
    end if

    ! If needed, evaluate (x1,x2) mixed derivatives at 4 corners
    if (self % bspline_2d % bs1 % bc_type == sll_p_hermite .and. &
        self % bspline_2d % bs2 % bc_type == sll_p_hermite) then
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

    ! Compute 2D spline that interpolates analytical 2D profile at points above
    call s_compute_interpolant_with_correct_arguments( &
      self % bspline_2d, &
      self % gtau, &
      val1_min    = bc1_min, &
      val1_max    = bc1_max, &
      val2_min    = bc2_min, &
      val2_max    = bc2_max, &
      val_corners = val_corners )

    ! Deallocate local arrays
    if (allocated(bc1_min))     deallocate( bc1_min )
    if (allocated(bc1_max))     deallocate( bc1_max )
    if (allocated(bc2_min))     deallocate( bc2_min )
    if (allocated(bc2_max))     deallocate( bc2_max )
    if (allocated(val_corners)) deallocate( val_corners )

  end subroutine init

  !-----------------------------------------------------------------------------
  subroutine run_tests( self, verbose, successful )

    class(t_bspline_2d_test_facility), intent(inout) :: self    ! TODO verify intent
    logical                          , intent(in   ) :: verbose
    logical                          , intent(  out) :: successful

    type(sll_t_bspline_2d ) :: bspline_2d
    type(t_profile_2d_info) :: info
    type(t_test_report    ) :: report
    integer                 :: i1, i2
    real(wp), pointer       :: tau1(:)
    real(wp), pointer       :: tau2(:)
    real(wp)                :: error
    real(wp)                :: max_norm_error

    successful = .false.

    ! Extract information about 2D analytical profile
    call self % profile_2d % get_info( info )

    ! Get spline interpolation points
    tau1 => self % bspline_2d % bs1 % tau
    tau2 => self % bspline_2d % bs2 % tau

    ! Print report to terminal on a single line
    write(*,'(6i10)', advance='no') &
      self%nx1, self%nx2, self%deg1, self%deg2, self%bc1, self%bc2

    ! Evaluate 2D spline at interpolation points: error should be zero
    max_norm_error = 0.0_wp
    do i2 = 1, size( tau2 )
      do i1 = 1, size( tau1 )
        error = self % gtau(i1,i2) &
              - sll_f_interpolate_value_2d( self % bspline_2d, tau1(i1), tau2(i2) )
        max_norm_error = max( max_norm_error, abs( error ) )
      end do
    end do

!    write(*,*) "Evaluate 2D spline at interpolation points: error should be zero"
!    write(*,*) "max_norm_error = ", max_norm_error
!    write(*,*)

    ! Store test data into report
    report % name   = 'zero_error'
    report % tol    = 1e-14_wp
    report % error  = max_norm_error
    report % passed = (report%error <= report%tol)

    ! Print report to terminal on a single line
    write(*,'(e12.2,L12)') report%error, report%passed

    ! Free spline memory
    ! TODO: maybe do it in another function
    call sll_s_bspline_2d_free( self % bspline_2d )

    ! Determine if all tests were successful
    successful = report % passed

  end subroutine run_tests

  !-----------------------------------------------------------------------------
  subroutine s_compute_interpolant_with_correct_arguments( &
      bspline_2d, &
      gtau      , &
      val1_min  , &
      val1_max  , &
      val2_min  , &
      val2_max  , &
      val_corners )

    type(sll_t_bspline_2d), intent(inout) :: bspline_2d
    real(wp)              , intent(in)    :: gtau(:,:)
    real(wp),     optional, intent(in)    :: val1_min(:,:)
    real(wp),     optional, intent(in)    :: val1_max(:,:)
    real(wp),     optional, intent(in)    :: val2_min(:,:)
    real(wp),     optional, intent(in)    :: val2_max(:,:)
    real(wp),     optional, intent(in)    :: val_corners(:,:,:)

    ! Hermite - other
    if (bspline_2d % bs1 % bc_type == sll_p_hermite .and. &
        bspline_2d % bs2 % bc_type /= sll_p_hermite) then

      call sll_s_compute_bspline_2d( bspline_2d, gtau, &
        val1_min = val1_min, &
        val1_max = val1_max )

    ! other - Hermite
    else if (bspline_2d % bs1 % bc_type /= sll_p_hermite .and. &
             bspline_2d % bs2 % bc_type == sll_p_hermite) then

      call sll_s_compute_bspline_2d( bspline_2d, gtau, &
        val2_min = val2_min, &
        val2_max = val2_max )

    ! Hermite - Hermite
    else if (bspline_2d % bs1 % bc_type == sll_p_hermite .and. &
             bspline_2d % bs2 % bc_type == sll_p_hermite) then

      call sll_s_compute_bspline_2d( bspline_2d, gtau, &
        val1_min    = val1_min, &
        val1_max    = val1_max, &
        val2_min    = val2_min, &
        val2_max    = val2_max, &
        val_corners = val_corners )

    ! other - other
    else

      call sll_s_compute_bspline_2d( bspline_2d, gtau )

    end if

  end subroutine s_compute_interpolant_with_correct_arguments

end module m_test_bsplines_2d
