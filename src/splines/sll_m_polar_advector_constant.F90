module sll_m_polar_advector_constant
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_polar_advector_base, only: sll_c_polar_advector

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_polar_advector) :: sll_t_polar_advector_constant

    ! Constant advection velocities
    real(wp) :: a(2)

  contains

    procedure :: init           => s_polar_advector_constant__init
    procedure :: velocity_field => s_polar_advector_constant__velocity_field
    procedure :: flow_field     => f_polar_advector_constant__flow_field
    procedure :: free           => s_polar_advector_constant__free

  end type sll_t_polar_advector_constant

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_polar_advector_constant__init( self, a )
    class(sll_t_polar_advector_constant), intent(inout) :: self
    real(wp)                            , intent(in   ) :: a(2) ! constant advection velocities

    self % a(:) = a(:)

  end subroutine s_polar_advector_constant__init

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_polar_advector_constant__velocity_field( self, x, a )
    class(sll_t_polar_advector_constant), intent(in   ) :: self
    real(wp)                            , intent(in   ) :: x(2)
    real(wp)                            , intent(  out) :: a(2)

    a(:) = self % a(:)

  end subroutine s_polar_advector_constant__velocity_field

  !-----------------------------------------------------------------------------
  SLL_PURE function f_polar_advector_constant__flow_field( self, x, h ) result( y )
    class(sll_t_polar_advector_constant), intent(in) :: self
    real(wp)                            , intent(in) :: x(2)
    real(wp)                            , intent(in) :: h
    real(wp) :: y(2)

    y(:) = x(:) + self % a(:) * h

  end function f_polar_advector_constant__flow_field

  !-----------------------------------------------------------------------------
  subroutine s_polar_advector_constant__free( self )
    class(sll_t_polar_advector_constant), intent(inout) :: self
  end subroutine s_polar_advector_constant__free

end module sll_m_polar_advector_constant
