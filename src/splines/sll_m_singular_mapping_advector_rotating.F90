module sll_m_singular_mapping_advector_rotating
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_singular_mapping_advector_base, only: sll_c_singular_mapping_advector

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type, extends(sll_c_singular_mapping_advector) :: sll_t_singular_mapping_advector_rotating

    real(wp) :: xc(2) ! center of rotation
    real(wp) :: omega ! angular velocity

  contains

    procedure :: init           => s_singular_mapping_advector_rotating__init
    procedure :: velocity_field => s_singular_mapping_advector_rotating__velocity_field
    procedure :: flow_field     => f_singular_mapping_advector_rotating__flow_field
    procedure :: free           => s_singular_mapping_advector_rotating__free

  end type sll_t_singular_mapping_advector_rotating

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_singular_mapping_advector_rotating__init( self, xc, omega )
    class(sll_t_singular_mapping_advector_rotating), intent(inout) :: self
    real(wp)                                       , intent(in   ) :: xc(2) ! center of rotation
    real(wp)                                       , intent(in   ) :: omega ! angular velocity

    self % xc(:) = xc(:)
    self % omega = omega

  end subroutine s_singular_mapping_advector_rotating__init

  !-----------------------------------------------------------------------------
  SLL_PURE subroutine s_singular_mapping_advector_rotating__velocity_field( self, x, a )
    class(sll_t_singular_mapping_advector_rotating), intent(in   ) :: self
    real(wp)                                       , intent(in   ) :: x(2)
    real(wp)                                       , intent(  out) :: a(2)

    a(1) = - self % omega * ( x(2) - self % xc(2) )
    a(2) =   self % omega * ( x(1) - self % xc(1) )

  end subroutine s_singular_mapping_advector_rotating__velocity_field

  !-----------------------------------------------------------------------------
  SLL_PURE function f_singular_mapping_advector_rotating__flow_field( self, x, h ) result( y )
    class(sll_t_singular_mapping_advector_rotating), intent(in) :: self
    real(wp)                                       , intent(in) :: x(2)
    real(wp)                                       , intent(in) :: h
    real(wp) :: y(2)

    y(1) =   self % xc(1) &
           + ( x(1) - self % xc(1) ) * cos( self % omega * h ) &
           - ( x(2) - self % xc(2) ) * sin( self % omega * h )

    y(2) =   self % xc(2) &
           + ( x(1) - self % xc(1) ) * sin( self % omega * h ) &
           + ( x(2) - self % xc(2) ) * cos( self % omega * h )

  end function f_singular_mapping_advector_rotating__flow_field

  !-----------------------------------------------------------------------------
  subroutine s_singular_mapping_advector_rotating__free( self )
    class(sll_t_singular_mapping_advector_rotating), intent(inout) :: self
  end subroutine s_singular_mapping_advector_rotating__free

end module sll_m_singular_mapping_advector_rotating
