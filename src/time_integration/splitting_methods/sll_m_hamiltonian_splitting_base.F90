module sll_m_hamiltonian_splitting_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_hamiltonian_splitting_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, abstract :: sll_c_hamiltonian_splitting_base

   contains 
     procedure(splitting), deferred :: lie_splitting
     procedure(splitting), deferred :: strang_splitting
     procedure(empty),     deferred :: free
  end type sll_c_hamiltonian_splitting_base


  abstract interface
     subroutine splitting(self, dt, number_steps)
       use sll_m_working_precision
       import sll_c_hamiltonian_splitting_base
       class(sll_c_hamiltonian_splitting_base), intent(inout)   :: self !< time splitting object
       sll_real64,                              intent(in)      :: dt !< time step size
       sll_int32,                               intent(in)      :: number_steps !< number of time steps
     end subroutine splitting
  end interface

  abstract interface
     subroutine empty(self)
       import sll_c_hamiltonian_splitting_base
       class(sll_c_hamiltonian_splitting_base), intent(inout)   :: self !< time splitting object
     end subroutine empty
  end interface

end module sll_m_hamiltonian_splitting_base
