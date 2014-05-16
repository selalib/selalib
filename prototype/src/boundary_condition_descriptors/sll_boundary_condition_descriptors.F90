!> @namespace sll_boundary_condition_descriptors
!> @brief Describe different boundary conditions
!> @details
!> The intent of this module is to provide a single, library-wide definition
!> of the names used to describe different boundary conditions. One should 
!> ALWAYS refer to specific boundary conditions by their
!> names and not through their integer representation, which could be changed.
!>
!> <b> How to use-it </b>
!>
!> Just add the line
!> @code
!> #include "sll_boundary_condition_descriptors.h"
!> @endcode
!
! To be considered here is to include also BC combinations, which may help
! save some coding instead of managing this internally within routines, for
! example a flag like SLL_DIRICHLET_NEUMANN could indicate two BC's along
! a particular dimension...
module sll_boundary_condition_descriptors
#include "sll_working_precision.h"
  implicit none

  sll_int32, parameter :: SLL_USER_DEFINED   = -1 
  sll_int32, parameter :: SLL_PERIODIC       = 0 
  sll_int32, parameter :: SLL_DIRICHLET      = 1 
  sll_int32, parameter :: SLL_NEUMANN        = 2
  sll_int32, parameter :: SLL_HERMITE        = 3
  sll_int32, parameter :: SLL_NEUMANN_MODE_0 = 4
  sll_int32, parameter :: SLL_SET_TO_LIMIT   = 5
  sll_int32, parameter :: SLL_INTERIOR       = 99

end module sll_boundary_condition_descriptors
