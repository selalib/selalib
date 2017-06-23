!> @authors Yaman Güçlü, IPP Garching
!>
!> @brief
!> f(x,y) = cos(2\pi x1) cos(2\pi x2)

module m_analytical_profiles_2d_cos_cos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: &
    f64

  use sll_m_constants, only: &
    sll_p_twopi

  use m_analytical_profiles_2d_base, only: &
    t_profile_2d_info, &
    c_analytical_profile_2d

  implicit none

  public :: &
    t_analytical_profile_2d_cos_cos

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Concrete type for 2D interpolation test-cases
  type, extends( c_analytical_profile_2d ) :: t_analytical_profile_2d_cos_cos

    ! Profile info (not overwritable by user)
    type(t_profile_2d_info), private :: info

  contains

    ! Constructor
    procedure :: init

    ! Abstract interface
    procedure :: get_info
    procedure :: eval

  end type t_analytical_profile_2d_cos_cos
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> Initialize profile
  !-----------------------------------------------------------------------------
  subroutine init( self )
    class( t_analytical_profile_2d_cos_cos ), intent(out) :: self

    ! Initialize basic profile info
    self % info % x1_min        = 0.0_wp
    self % info % x1_max        = 1.0_wp
    self % info % x2_min        = 0.0_wp
    self % info % x2_max        = 1.0_wp
    self % info % x1_periodic   = .true.
    self % info % x2_periodic   = .true.
    self % info % x1_poly_order = -1
    self % info % x2_poly_order = -1

  end subroutine init

  !-----------------------------------------------------------------------------
  !> Get profile info
  !-----------------------------------------------------------------------------
  pure subroutine get_info( self, info )
    class( t_analytical_profile_2d_cos_cos ), intent(in   ) :: self
    type ( t_profile_2d_info               ), intent(  out) :: info

    info = self % info

  end subroutine get_info

  !-----------------------------------------------------------------------------
  !> Evaluate 2D profile (or one of its derivatives)
  !-----------------------------------------------------------------------------
  pure function eval( self, x1, x2, diff_x1, diff_x2 ) result( f )
    class( t_analytical_profile_2d_cos_cos ), intent(in) :: self
    real(wp)                                , intent(in) :: x1
    real(wp)                                , intent(in) :: x2
    integer,                        optional, intent(in) :: diff_x1
    integer,                        optional, intent(in) :: diff_x2
    real(wp) :: f

    integer  :: d1  , d2
    real(wp) :: arg1, arg2
    real(wp) :: f1  , f2

    if (present( diff_x1 )) then; d1 = diff_x1; else; d1 = 0; end if
    if (present( diff_x2 )) then; d2 = diff_x2; else; d2 = 0; end if

    arg1 = sll_p_twopi * x1
    arg2 = sll_p_twopi * x2

    select case (modulo( d1, 4 ))
    case( 0 )
      f1 =  cos( arg1 )
    case( 1 )
      f1 = -sin( arg1 )
    case( 2 )
      f1 = -cos( arg1 )
    case( 3 )
      f1 =  sin( arg1 )
    end select

    select case (modulo( d2, 4 ))
    case( 0 )
      f2 =  cos( arg2 )
    case( 1 )
      f2 = -sin( arg2 )
    case( 2 )
      f2 = -cos( arg2 )
    case( 3 )
      f2 =  sin( arg2 )
    end select

    f = sll_p_twopi**(d1+d2) * f1 * f2

  end function eval



end module m_analytical_profiles_2d_cos_cos
