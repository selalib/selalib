!> @ingroup time_integration
!> @author Katharina Kormann, IPP
!> @brief Base class for Hamiltonian splittings.
!> @details  Contains deferred function for strang splitting, lie splitting and lie splitting with oposite ordering of the split-steps.
!> Moreover, composition methods based on lie and strang splitting are implemented (cf. Hairer, Lubich, Wanner, Geometric numeric integration.
module sll_m_time_propagator_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

use sll_m_hamiltonian_splitting_base


  implicit none

  public :: &
    sll_c_time_propagator_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Type for Hamiltonian splittings
  type, abstract, extends(sll_c_hamiltonian_splitting_base) :: sll_c_time_propagator_base

  
  end type sll_c_time_propagator_base


  
end module sll_m_time_propagator_base
