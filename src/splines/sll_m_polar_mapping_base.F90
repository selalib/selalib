module sll_m_polar_mapping_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  implicit none

  public :: sll_c_polar_mapping

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract  type, polar mapping
  type, abstract :: sll_c_polar_mapping

  contains

    ! Deferred procedures
    procedure(i_fun_eval)      , deferred :: eval
    procedure(i_fun_jmat)      , deferred :: jmat ! Jacobian matrix
    procedure(i_sub_store_data), deferred :: store_data

    ! Non-deferred procedures
    procedure :: jdet => f_polar_mapping__jdet ! Jacobian determinant

  end type sll_c_polar_mapping

  ! Interfaces for deferred procedures
  abstract interface

    SLL_PURE function i_fun_eval( self, eta ) result( x )
      import sll_c_polar_mapping, wp
      class(sll_c_polar_mapping), intent(in) :: self
      real(wp)                  , intent(in) :: eta(2)
      real(wp) :: x(2)
    end function i_fun_eval

    SLL_PURE function i_fun_jmat( self, eta ) result( jmat )
      import sll_c_polar_mapping, wp
      class(sll_c_polar_mapping), intent(in) :: self
      real(wp)                  , intent(in) :: eta(2)
      real(wp) :: jmat(2,2)
    end function i_fun_jmat

    subroutine i_sub_store_data( self, n1, n2, file_name )
      import sll_c_polar_mapping
      class(sll_c_polar_mapping), intent(in) :: self
      integer                   , intent(in) :: n1
      integer                   , intent(in) :: n2
      character(len=*)          , intent(in) :: file_name
    end subroutine i_sub_store_data

  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SLL_PURE function f_polar_mapping__jdet( self, eta ) result( jdet )
    class(sll_c_polar_mapping), intent(in) :: self
    real(wp)                  , intent(in) :: eta(2)
    real(wp) :: jdet

    real(wp) :: jmat(2,2)

    jmat = self % jmat( eta )

    jdet = jmat(1,1) * jmat(2,2) - jmat(1,2) * jmat(2,1)

  end function f_polar_mapping__jdet

end module sll_m_polar_mapping_base
