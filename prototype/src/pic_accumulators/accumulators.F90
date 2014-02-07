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

  type field_accumulator_cell! for interpolating the elec field from
                     ! its grid values(in this type) in the particles
     sll_real64 :: Ex_sw
     sll_real64 :: Ex_se
     sll_real64 :: Ex_nw
     sll_real64 :: Ex_ne
     sll_real64 :: Ey_sw
     sll_real64 :: Ey_se
     sll_real64 :: Ey_nw
     sll_real64 :: Ey_ne
  end type field_accumulator_cell


contains
  
  subroutine sll_calculate_chargedensity( &
                cells_number,  &
                p_group, &
                charge_density )

    sll_int32, intent(in) :: cells_number
    type(sll_particle_group_2d), intent(in) :: p_group
    type(charge_accumulator_cell), dimension(:), intent(out) :: charge_density
    sll_int64    ::  j
    
    charge_density = 0._f64
    do j = 1, p_group%num_particles
       charge_density(p_group%p_list(j)%ic)%q_sw = &
            charge_density(p_group%p_list(j)%ic)%q_sw + &
            p_group%p_list(j)%q &
            * (1._f64 - p_group%p_list(j)%dx/m2d%delta_eta1) &
            * (1._f64 - p_group%p_list(j)%dy/m2d%delta_eta2)
       charge_density(p_group%p_list(j)%ic)%q_se = &

       charge_density(p_group%p_list(j)%ic)%q_nw = &

       charge_density(p_group%p_list(j)%ic)%q_ne = &
            
    enddo

  end subroutine sll_calculate_chargedensity
  
end module sll_accumulators
