!> @authors  Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> f(x,y) is polynomial in x and y

module m_analytical_profiles_2d_poly
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

  use sll_m_working_precision, only: &
    f64

  use sll_m_constants, only: &
    sll_p_twopi

  use m_analytical_profiles_2d_base, only: &
    t_profile_2d_info, &
    c_analytical_profile_2d

  implicit none

  public :: &
    t_analytical_profile_2d_poly

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Concrete type for 2D interpolation test-cases
  type, extends( c_analytical_profile_2d ) :: t_analytical_profile_2d_poly

    ! Profile info (not overwritable by user)
    type(t_profile_2d_info), private :: info

    integer :: deg1
    integer :: deg2
    real(wp), allocatable :: coeffs(:,:)

  contains

    ! Constructor
    procedure :: init

    ! Abstract interface
    procedure :: get_info
    procedure :: eval
    procedure :: max_norm ! not implemented for this profile

  end type t_analytical_profile_2d_poly
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !-----------------------------------------------------------------------------
  !> Initialize profile
  !-----------------------------------------------------------------------------
  subroutine init( self, deg1, deg2 )
    class( t_analytical_profile_2d_poly ), intent(  out) :: self
    integer                              , intent(in   ) :: deg1
    integer                              , intent(in   ) :: deg2

    ! Initialize basic profile info
    self % info % x1_min        = -1.0_wp
    self % info % x1_max        =  1.0_wp
    self % info % x2_min        = -1.0_wp
    self % info % x2_max        =  1.0_wp
    self % info % x1_periodic   = .false.
    self % info % x2_periodic   = .false.
    self % info % x1_poly_order =  deg1
    self % info % x2_poly_order =  deg2

    ! Generate random polynomial coefficients
    if( allocated( self%coeffs ) ) deallocate( self%coeffs )
    allocate( self % coeffs (0:deg1,0:deg2) )
    self % deg1 = deg1
    self % deg2 = deg2
    call random_number( self % coeffs )     ! 0 <= c < 1
    self % coeffs = 1.0_wp - self % coeffs  ! 0 < c <= 1
    self % coeffs = 0.2_wp + 0.8_wp * self % coeffs ! 0.2 < c <= 1

  end subroutine init

  !-----------------------------------------------------------------------------
  !> Get profile info
  !-----------------------------------------------------------------------------
  pure subroutine get_info( self, info )
    class( t_analytical_profile_2d_poly ), intent(in   ) :: self
    type ( t_profile_2d_info            ), intent(  out) :: info

    info = self % info

  end subroutine get_info

  !-----------------------------------------------------------------------------
  !> Evaluate 2D profile (or one of its derivatives)
  !-----------------------------------------------------------------------------
  pure function eval( self, x1, x2, diff_x1, diff_x2 ) result( f )
    class( t_analytical_profile_2d_poly ), intent(in) :: self
    real(wp)                             , intent(in) :: x1
    real(wp)                             , intent(in) :: x2
    integer,                     optional, intent(in) :: diff_x1
    integer,                     optional, intent(in) :: diff_x2
    real(wp) :: f

    integer  :: d1, d2
    integer  :: i1, i2
    real(wp) :: c1, c2

    if (present(diff_x1)) then; d1 = diff_x1; else; d1 = 0; end if
    if (present(diff_x2)) then; d2 = diff_x2; else; d2 = 0; end if

    f = 0.0_wp
    if (d1 > self%deg1 .or. d2 > self%deg2) return

    do i1 = d1, self%deg1
      c1 = real( falling_factorial( i1, d1 ), kind=wp ) * x1**(i1-d1)
      do i2 = d2, self%deg2
        c2 = real( falling_factorial( i2, d2 ), kind=wp ) * x2**(i2-d2)
        f  = f + self%coeffs(i1,i2) * c1 * c2
      end do
    end do

  end function eval

  !-----------------------------------------------------------------------------
  !> Evaluate max norm of profile (or one of its derivatives) over domain
  !-----------------------------------------------------------------------------
  function max_norm( self, diff_x1, diff_x2 ) result( norm )
    class( t_analytical_profile_2d_poly ), intent(in) :: self
    integer,                     optional, intent(in) :: diff_x1
    integer,                     optional, intent(in) :: diff_x2
    real(wp) :: norm

    integer  :: d1, d2

    if (present(diff_x1)) then; d1 = diff_x1; else; d1 = 0; end if
    if (present(diff_x2)) then; d2 = diff_x2; else; d2 = 0; end if

    if (self%info%x1_max >= abs( self%info%x1_min ) .and. &
        self%info%x2_max >= abs( self%info%x2_min )) then
      ! For x1_max >= |x1_min| and x2_max >= |x2_min|:
      ! max(|f^(d1,d2)(x1,x2)|) = f^(d1,d2)(x1_max,x2_max)
      norm = self % eval( self%info%x1_max, self%info%x2_max, d1, d2 )
    else
      SLL_ERROR( "t_analytical_profile_2d_poly % max_norm", "General formula not implemented" )
    end if

  end function max_norm

  !-----------------------------------------------------------------------------
  !> Calculate falling factorial of x
  !> [ https://en.wikipedia.org/wiki/Falling_and_rising_factorials ]
  !-----------------------------------------------------------------------------
  pure function falling_factorial( x, n ) result( c )
    integer, intent(in) :: x
    integer, intent(in) :: n
    integer :: c
    integer :: k
    c = 1
    do k = 0, n-1
      c = c * (x-k)
    end do
  end function falling_factorial

end module m_analytical_profiles_2d_poly
