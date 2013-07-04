module sll_boundary_condition_descriptors
#include "sll_working_precision.h"
  implicit none

  ! The intent of this module is to provide a single, library-wide definition
  ! of the names used to describe different boundary conditions. One should 
  ! ALWAYS refer to specific boundary conditions by their
  ! names and not through their integer representation, which could be changed.
  !
  ! To be considered here is to include also BC combinations, which may help
  ! save some coding instead of managing this internally within routines, for
  ! example a flag like SLL_DIRICHLET_NEUMANN could indicate two BC's along
  ! a particular dimension...

  sll_int32, parameter :: SLL_PERIODIC  = 0
  sll_int32, parameter :: SLL_DIRICHLET = 1 
  sll_int32, parameter :: SLL_NEUMANN   = 2
  sll_int32, parameter :: SLL_HERMITE   = 3

end module sll_boundary_condition_descriptors
