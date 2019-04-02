module m_splines_error_bounds
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision , only: f64
  use m_analytical_profiles_1d, only: c_analytical_profile_1d
  use m_analytical_profiles_2d, only: c_analytical_profile_2d

  public :: &
    sll_f_spline_1d_error_bound         , &
    sll_f_spline_1d_error_bound_on_deriv, &
    sll_f_spline_2d_error_bound         , &
    sll_f_spline_2d_error_bounds_on_grad, &
    sll_f_spline_2d_error_bound_on_mixed_deriv

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  !
  !> @param[in] h       cell width
  !> @param[in] deg     degree of spline S
  !> @param[in] norm_f  max(f) (or its derivative) over domain
  !> @result    norm_e  max of error $E(x):=f(x)-S(x)$ over domain
  !-----------------------------------------------------------------------------
  pure function f_tihomirov_error_bound( h, deg, norm_f ) result( norm_e )
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

  end function f_tihomirov_error_bound

  !-----------------------------------------------------------------------------
  function sll_f_spline_1d_error_bound( profile_1d, dx, deg ) result( max_error )
    class(c_analytical_profile_1d), intent(in) :: profile_1d
    real(wp)                      , intent(in) :: dx
    integer                       , intent(in) :: deg
    real(wp) :: max_error

    real(wp) :: max_norm
    max_norm  = profile_1d % max_norm( deg+1 )
    max_error = f_tihomirov_error_bound( dx, deg, max_norm )

  end function sll_f_spline_1d_error_bound

  !-----------------------------------------------------------------------------
  function sll_f_spline_1d_error_bound_on_deriv( profile_1d, dx, deg ) result( max_error )
    class(c_analytical_profile_1d), intent(in) :: profile_1d
    real(wp)                      , intent(in) :: dx
    integer                       , intent(in) :: deg
    real(wp) :: max_error

    real(wp) :: max_norm
    max_norm  = profile_1d % max_norm( deg+1 )
    max_error = f_tihomirov_error_bound( dx, deg-1, max_norm )

  end function sll_f_spline_1d_error_bound_on_deriv

  !-----------------------------------------------------------------------------
  function sll_f_spline_2d_error_bound( profile_2d, dx1, dx2, deg1, deg2 ) result( max_error )
    class(c_analytical_profile_2d), intent(in) :: profile_2d
    real(wp)                      , intent(in) :: dx1
    real(wp)                      , intent(in) :: dx2
    integer                       , intent(in) :: deg1
    integer                       , intent(in) :: deg2
    real(wp) :: max_error

    real(wp) :: max_norm1
    real(wp) :: max_norm2

    ! Max norm of highest partial derivatives in x1 and x2 of analytical profile
    max_norm1 = profile_2d % max_norm( deg1+1, 0      )
    max_norm2 = profile_2d % max_norm( 0     , deg2+1 )

    ! Error bound on function value
    max_error = f_tihomirov_error_bound( dx1, deg1, max_norm1 ) &
              + f_tihomirov_error_bound( dx2, deg2, max_norm2 )

    ! Empirical correction: for linear interpolation increase estimate by 5%
    if (deg1 == 1 .or. deg2 == 1) then
      max_error = 1.05_wp * max_error
    end if

  end function sll_f_spline_2d_error_bound

  !-----------------------------------------------------------------------------
  function sll_f_spline_2d_error_bounds_on_grad( profile_2d, dx1, dx2, deg1, deg2 ) result( max_error )
    class(c_analytical_profile_2d), intent(in) :: profile_2d
    real(wp)                      , intent(in) :: dx1
    real(wp)                      , intent(in) :: dx2
    integer                       , intent(in) :: deg1
    integer                       , intent(in) :: deg2
    real(wp) :: max_error(2)

    real(wp) :: max_norm1
    real(wp) :: max_norm2

    ! Max norm of highest partial derivatives in x1 and x2 of analytical profile
    max_norm1 = profile_2d % max_norm( deg1+1, 0      )
    max_norm2 = profile_2d % max_norm( 0     , deg2+1 )

    ! Error bound on x1-derivative
    max_error(1) = f_tihomirov_error_bound( dx1, deg1-1, max_norm1 ) &
                 + f_tihomirov_error_bound( dx2, deg2  , max_norm2 )

    ! Error bound on x2-derivative
    max_error(2) = f_tihomirov_error_bound( dx1, deg1  , max_norm1 ) &
                 + f_tihomirov_error_bound( dx2, deg2-1, max_norm2 )

  end function sll_f_spline_2d_error_bounds_on_grad

  !-----------------------------------------------------------------------------
  ! NOTE: The following estimate has no theoretical justification but captures
  !       the correct asympthotic rate of convergence.
  !       The error constant is probably overestimated.
  !-----------------------------------------------------------------------------
  function sll_f_spline_2d_error_bound_on_mixed_deriv( profile_2d, dx1, dx2, deg1, deg2 ) result( max_error )
    class(c_analytical_profile_2d), intent(in) :: profile_2d
    real(wp)                      , intent(in) :: dx1
    real(wp)                      , intent(in) :: dx2
    integer                       , intent(in) :: deg1
    integer                       , intent(in) :: deg2
    real(wp) :: max_error

    real(wp) :: max_norm1
    real(wp) :: max_norm2

    ! Max norm of highest partial derivatives in x1 and x2 of analytical profile
    max_norm1 = profile_2d % max_norm( deg1+1, 1      )
    max_norm2 = profile_2d % max_norm( 1     , deg2+1 )

    ! Error bound on x1-x2 mixed derivative
    max_error = f_tihomirov_error_bound( dx1, deg1-1, max_norm1 ) &
              + f_tihomirov_error_bound( dx2, deg2-1, max_norm2 )

  end function sll_f_spline_2d_error_bound_on_mixed_deriv

end module m_splines_error_bounds
