!> @ingroup boundary_condition_descriptors
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
!> #include "sll_m_boundary_condition_descriptors.h"
!> @endcode
!
! To be considered here is to include also BC combinations, which may help
! save some coding instead of managing this internally within routines, for
! example a flag like SLL_DIRICHLET_NEUMANN could indicate two BC's along
! a particular dimension...
module sll_m_boundary_condition_descriptors
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_conductor, &
    sll_dirichlet, &
    sll_d_halo, &
    sll_d_one_sided, &
    sll_hermite, &
    sll_interior, &
    sll_neumann, &
    sll_neumann_mode_0, &
    sll_periodic, &
    sll_set_to_limit, &
    sll_silver_muller, &
    sll_user_defined

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> User defined boundary condition
  sll_int32, parameter :: SLL_USER_DEFINED   = -1 
  !> Periodic boundary condition u(1)=u(n)
  sll_int32, parameter :: SLL_PERIODIC       = 0 
  !> Dirichlet boundary condition 
  sll_int32, parameter :: SLL_DIRICHLET      = 1 
  !> Neumann boundary condition 
  sll_int32, parameter :: SLL_NEUMANN        = 2
  !> Hermite boundary condition
  sll_int32, parameter :: SLL_HERMITE        = 3
  !> Neumann boundary condition
  sll_int32, parameter :: SLL_NEUMANN_MODE_0 = 4
  !> PLEASE ADD DOCUMENTATION
  sll_int32, parameter :: SLL_SET_TO_LIMIT   = 5
  !> Interior of domain
  sll_int32, parameter :: SLL_INTERIOR       = 6
  !> Incoming wave boundar condition for Maxwell
  sll_int32, parameter :: SLL_INCOMING_WAVE  = 7
  !> Metallic boundary condition for Maxwell
  sll_int32, parameter :: SLL_CONDUCTOR      = 8
  !> Absorbing boundary condition fro Maxwell
  sll_int32, parameter :: SLL_SILVER_MULLER  = 9
  !> Use a one-sided stencil at the boundary
  sll_int32, parameter :: SLL_D_ONE_SIDED    = 10
  !> Values outside the domain are provided as halo cells (for domain decomposition)
  sll_int32, parameter :: SLL_D_HALO         = 11


end module sll_m_boundary_condition_descriptors
