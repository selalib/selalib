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

module sll_lt_pic_particle_module
!module sll_lt_particle_representations  ! old name

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_cartesian_meshes

  implicit none

!> @typedef sll_lt_pic_4d_particle
!! particles of type sll_lt_pic_4d_particle know their neighbors in the cartesian grid used to initialize or remap the particles
!! ** in every non-periodic dimension: the first and last particles have no neighbors and by convention the ngb_index is the
!! particle own index.
!! ** in the periodic dimensions: the first and last particles are neighbors to each other, so that we can travel through the
!! 'periodic boundaries'

  type :: sll_lt_pic_4d_particle
!  type :: sll_lt_particle_4d  ! old name
     sll_int32  :: ic   ! cell index, linearly arranged
     sll_real32 :: dx
     sll_real32 :: dy
     sll_real64 :: vx
     sll_real64 :: vy
     sll_real32 :: q
     ! <<neighbour_pointers>> particle indices of logical neighbors
     sll_int64  :: ngb_xleft_index   ! 32 devrait suffire
     sll_int64  :: ngb_xright_index
     sll_int64  :: ngb_yleft_index
     sll_int64  :: ngb_yright_index
     sll_int64  :: ngb_vxleft_index
     sll_int64  :: ngb_vxright_index
     sll_int64  :: ngb_vyleft_index
     sll_int64  :: ngb_vyright_index               
  end type sll_lt_pic_4d_particle


  type :: sll_lt_pic_4d_particle_guard
     type(sll_lt_pic_4d_particle), pointer :: p
  end type sll_lt_pic_4d_particle_guard

  type :: sll_lt_pic_4d_particle_guard_ptr
     type(sll_lt_pic_4d_particle_guard), dimension(:), pointer :: g_list
  end type sll_lt_pic_4d_particle_guard_ptr


end module sll_lt_pic_particle_module
