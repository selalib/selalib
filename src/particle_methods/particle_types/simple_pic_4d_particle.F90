!**************************************************************
!  Copyright INRIA
!  Authors : 
!     ????
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

!> @ingroup particle_methods

!> @author MCP ALH

!> @brief Most simple particle model

!> @details Mainly useful to compare against more advanced models. Defined in
!> [[selalib:src/particle_methods/particle_types/simple_pic_4d_particle.F90]]. A group of sll_simple_pic_4d_particle is
!> a sll_simple_pic_4d_group
!> [[selalib:src/particle_methods/particle_types/simple_pic_4d_group.F90::sll_simple_pic_4d_group]]

module sll_simple_pic_4d_particle_module

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_cartesian_meshes

  implicit none

  type :: sll_simple_pic_4d_particle
     !> cell index in the x dimension (can be outside physical domain)
     sll_int32  :: i_cell_x

     !> cell index in the y dimension (can be outside physical domain)
     sll_int32  :: i_cell_y

     !> relative position in cell (in [0,1])
     sll_real32 :: offset_x

     !> relative position in cell (in [0,1])
     sll_real32 :: offset_y

     !> speed along X axis
     sll_real64 :: v_x

     !> speed along Y axis
     sll_real64 :: v_y

     sll_real32 :: weight

  end type sll_simple_pic_4d_particle

end module sll_simple_pic_4d_particle_module
