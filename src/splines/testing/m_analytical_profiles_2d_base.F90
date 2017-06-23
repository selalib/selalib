!> @authors Yaman Güçlü, IPP Garching
!>
!> @brief
!> Define 2D functions on logical rectangular domain to test interpolation
!>

module m_analytical_profiles_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  implicit none

  public :: &
    t_profile_2d_info, &
    c_analytical_profile_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !-----------------------------------------------------------------------------
  !> Structure with basic info about 2D profile
  !-----------------------------------------------------------------------------
  type :: t_profile_2d_info

    real(wp) :: x1_min, x1_max ! x1 limits
    real(wp) :: x2_min, x2_max ! x2 limits
    logical  :: x1_periodic    ! .true. if profile is periodic along x1
    logical  :: x2_periodic    ! .true. if profile is periodic along x2
    integer  :: x1_poly_order  ! polynomial order in x1 (-1 if not a polynomial)
    integer  :: x2_poly_order  ! polynomial order in x2 (-1 if not a polynomial)

  end type t_profile_2d_info

  !-----------------------------------------------------------------------------
  !> Abstract type for 2D analytical profile to be used for testing interpolation
  !-----------------------------------------------------------------------------
  type, abstract :: c_analytical_profile_2d

  contains
    ! Get profile info
    procedure( i_subr_get_info     ), deferred :: get_info

    ! Evaluate 2D profile (or one of its derivatives)
    procedure( i_func_eval_profile ), deferred :: eval

  end type c_analytical_profile_2d
 
  !-----------------------------------------------------------------------------
  abstract interface

    ! Get profile info
    pure subroutine i_subr_get_info( self, info )
      import c_analytical_profile_2d, t_profile_2d_info
      class( c_analytical_profile_2d ), intent(in   ) :: self
      type ( t_profile_2d_info       ), intent(  out) :: info
    end subroutine i_subr_get_info

    ! Evaluate 2D profile (or one of its derivatives)
    pure function i_func_eval_profile( self, x1, x2, diff_x1, diff_x2 ) result( f )
      import c_analytical_profile_2d, wp
      class( c_analytical_profile_2d ), intent(in) :: self
      real(wp)                        , intent(in) :: x1
      real(wp)                        , intent(in) :: x2
      integer,                optional, intent(in) :: diff_x1
      integer,                optional, intent(in) :: diff_x2
      real(wp) :: f
    end function i_func_eval_profile

  end interface

end module m_analytical_profiles_2d_base
