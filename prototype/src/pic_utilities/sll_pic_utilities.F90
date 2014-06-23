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

module sll_pic_utilities
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_accumulators.h" 
  use sll_particle_group_2d_module

contains

  subroutine sll_first_charge_accumulation_2d( p_group, q_accum )
    type(sll_particle_group_2d), pointer         :: p_group
    type(sll_charge_accumulator_2d), pointer     :: q_accum
    type(sll_particle_2d), dimension(:), pointer :: p
    sll_int64  :: i
    sll_int32  :: num_cells
    sll_int64  :: num_particles
    sll_real64 :: tmp1
    sll_real64 :: tmp2

    SLL_ASSERT( associated(p_group) .and. associated(q_accum))
    num_particles =  p_group%number_particles
    p             => p_group%p_list

    do i=1,num_particles
       SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp1,tmp2)
    end do
  end subroutine sll_first_charge_accumulation_2d

  subroutine sll_first_charge_accumulation_2d_CS( p_group, q_accum )
    type(sll_particle_group_2d), pointer         :: p_group
    type(sll_charge_accumulator_2d_CS), pointer     :: q_accum
    type(sll_particle_2d), dimension(:), pointer :: p
    sll_int64  :: i
    sll_int32  :: num_cells
    sll_int64  :: num_particles
    sll_real64 :: tmp(1:4,1:2), temp

    SLL_ASSERT( associated(p_group) .and. associated(q_accum))
    num_particles =  p_group%number_particles
    p             => p_group%p_list

    do i=1,num_particles
       SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum,p(i),tmp,temp)
    end do
  end subroutine sll_first_charge_accumulation_2d_CS


end module sll_pic_utilities
