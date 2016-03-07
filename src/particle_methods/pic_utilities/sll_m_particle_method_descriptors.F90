 !> @ingroup pic_utilities
 !> @author Jakob Ameres
 !> @brief Descriptors for particle methods
 !> @details Should be replaced by sll_m_descriptors
module sll_m_particle_method_descriptors
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_p_collisions_none, &
    sll_p_controlvariate_maxwellian, &
    sll_p_controlvariate_none, &
    sll_p_controlvariate_standard, &
    sll_p_moment_match_initial, &
    sll_p_moment_match_none

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !-------------- CONTROL VARIATE / DELTA F -------------------------
  !> Descriptors concerning simulations with control variate
  
  !> No Control Variate - also known as full-f
  sll_int32, parameter :: sll_p_controlvariate_none = 0
  
  !> standard control variate for simulation
  sll_int32, parameter :: sll_p_controlvariate_standard=1
  
  !>Maxwellian control variate
  sll_int32, parameter :: sll_p_controlvariate_maxwellian=2
  
  !>Local Maxwellian control variate, estimates temperature etc.
  sll_int32, parameter :: SLL_CONTROLVARIATE_MAXWELLIAN_LOCAL=3
  
  !>use initial distribution f(t=0,x,v) as control variate
  sll_int32, parameter :: SLL_CONTROLVARIATE_INITIAL=4

!  character(len=*), parameter :: &
!          sll_controlvariate_key(0:4) = &
!         (/"sll_p_controlvariate_none            ",&
!           "sll_p_controlvariate_standard        ",&
!           "sll_p_controlvariate_maxwellian      ",&
!           "SLL_CONTROLVARIATE_MAXWELLIAN_LOCAL",&
!           "SLL_CONTROLVARIATE_INITIAL         "/)
!  
  
  
  !-------------- RANDOM NUMBERS -----------------------------
  !> No Control Variate - also known as full-f
  sll_int32, parameter :: SLL_HAMMERSLEY = 1
  
  sll_int32, parameter :: SLL_RANDOM_SOBOL = 2
  
  sll_int32, parameter :: SLL_RANDOM_SOBOL_SCRAMBLED = 3
  
  sll_int32, parameter :: SLL_RANDOM_HAMMERSLEY = 4
 
 
  !---------------- MOMENT MATCHING ---------------------------
  sll_int32, parameter :: sll_p_moment_match_none = 0
  !>Match initial moments
  sll_int32, parameter :: sll_p_moment_match_initial = 1
  !>match the moments of the sampling density
  sll_int32, parameter :: SLL_MOMENT_MATCH_PRIOR = 2
  !>match the moments of the sampling density, only velocity
  sll_int32, parameter :: SLL_MOMENT_MATCH_PRIOR_V = 3
  !>match the moments of the sampling density, only spatial
  sll_int32, parameter :: SLL_MOMENT_MATCH_PRIOR_X = 4
 
 
 
  
  !------------- COLLISION OPERATORS -------------------------
  !> do not implement collisions
  sll_int32, parameter :: sll_p_collisions_none = 0
  
  !>use the standard operator
  sll_int32, parameter :: SLL_COLLISIONS_STANDARD = 1
  
  !>use a krook operator for collisions
  sll_int32, parameter :: SLL_COLLISIONS_KROOK = 2
  
  !>use a landau operator for collisions
  sll_int32, parameter :: SLL_COLLISIONS_LANDAU = 3
  
  !>Key for collision operators
!   character(len=*), parameter :: &
!          sll_collisions_key(0:3) = &
!         (/"sll_p_collisions_none       ",&
!           "SLL_COLLISIONS_STANDARD   ",&
!           "SLL_COLLISIONS_KROOK      ",&
!           "SLL_COLLISIONS_LANDAU     " /) 
  
end module sll_m_particle_method_descriptors
 
