!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************


module sll_accumulators
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_particle_representations
  implicit none
  
  type charge_accumulator_cell! for particles deposition on the grid
     sll_real64 :: q_sw
     sll_real64 :: q_se
     sll_real64 :: q_nw
     sll_real64 :: q_ne
  end type charge_accumulator_cell

  type field_accumulator_cell! for interpolating theelec field from
     sll_real64 :: Ex_sw
     sll_real64 :: Ex_se
     sll_real64 :: Ex_nw
     sll_real64 :: Ex_ne
     sll_real64 :: Ey_sw
     sll_real64 :: Ey_se
     sll_real64 :: Ey_nw
     sll_real64 :: Ey_ne
                     ! its grid values(in this type) in the paticles
  end type field_accumulator_cell


contains
  
  subroutine calculate_chargedensity( &
                cells_number,  &
                particle_list, &
                charge_density )

    sll_int32, intent(in) :: cells_number
    type(sll_particle_group_2d), pointer, intent(in) :: particle_list
    type(charge_accumulator_cell), dimension(:,:), intent(out) :: charge_density

  end subroutine calculate_chargedensity
  
end module sll_accumulators
