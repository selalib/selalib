module sll_m_linear_operator_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_vector_space_base, only: sll_c_vector_space

  implicit none

  public :: sll_c_linear_operator

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: sll_c_linear_operator

  contains

    ! Deferred procedures
    procedure(i_fun_shape), deferred :: get_shape
    procedure(i_sub_dot ) , deferred :: dot
    procedure(i_sub_free) , deferred :: free

  end type sll_c_linear_operator

  ! Interfaces for deferred procedures
  abstract interface

    function i_fun_shape( self ) result( s )
      import sll_c_linear_operator
      class(sll_c_linear_operator), intent(in) :: self
      integer :: s(2)
    end function i_fun_shape

    subroutine i_sub_dot( self, x, y )
      import sll_c_linear_operator, sll_c_vector_space
      class(sll_c_linear_operator), intent(in   ) :: self
      class(sll_c_vector_space)   , intent(in   ) :: x
      class(sll_c_vector_space)   , intent(inout) :: y
    end subroutine i_sub_dot

    subroutine i_sub_free( self )
      import sll_c_linear_operator
      class(sll_c_linear_operator), intent(inout) :: self
    end subroutine i_sub_free

  end interface

end module sll_m_linear_operator_base
