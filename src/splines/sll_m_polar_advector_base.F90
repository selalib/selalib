module sll_m_polar_advector_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, abstract :: sll_c_polar_advector

  contains

    ! Deferred procedures
    procedure(i_sub_velocity_field), deferred :: velocity_field ! RHS of characteristic equations
    procedure(i_fun_flow_field)    , deferred :: flow_field ! analytical characteristics
    procedure(i_sub_free)          , deferred :: free

    ! Non-deferred procedures
    procedure :: advect => f_polar_advector__advect ! time integrator

  end type sll_c_polar_advector

  ! Interfaces for deferred procedures
  abstract interface

    subroutine i_sub_velocity_field( self, x, a )
      import sll_c_polar_advector, wp
      class(sll_c_polar_advector), intent(in   ) :: self
      real(wp)                   , intent(in   ) :: x(2)
      real(wp)                   , intent(  out) :: a(2)
    end subroutine i_sub_velocity_field

    SLL_PURE function i_fun_flow_field( self, x, h ) result( y )
      import sll_c_polar_advector, wp
      class(sll_c_polar_advector), intent(in) :: self
      real(wp)                   , intent(in) :: x(2)
      real(wp)                   , intent(in) :: h
      real(wp) :: y(2)
    end function i_fun_flow_field

    subroutine i_sub_free( self )
      import sll_c_polar_advector
      class(sll_c_polar_advector), intent(inout) :: self
    end subroutine i_sub_free

  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Explicit Runge-Kutta 4th order
  SLL_PURE function f_polar_advector__advect( self, x, h ) result( y )
    class(sll_c_polar_advector), intent(in) :: self
    real(wp)                   , intent(in) :: x(2)
    real(wp)                   , intent(in) :: h
    real(wp) :: y(2)

    real(wp) :: k1(2), k2(2), k3(2), k4(2)

    call self % velocity_field( x                  , k1 )
    call self % velocity_field( x + 0.5_wp * h * k1, k2 )
    call self % velocity_field( x + 0.5_wp * h * k2, k3 )
    call self % velocity_field( x +          h * k3, k4 )

    y = x + h * ( k1 + 2.0_wp * k2 + 2.0_wp * k3 + k4 ) / 6.0_wp

  end function f_polar_advector__advect

end module sll_m_polar_advector_base
