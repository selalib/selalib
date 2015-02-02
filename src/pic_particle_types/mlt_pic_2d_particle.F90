!**************************************************************
!  Copyright INRIA
!  Authors : MCP,ALH
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

!> @file
!> @brief Data structures representing one multilevel linearly transformed particle

! <<mlt_pic_naming_convention>> The following name parts are sorted from the most global to the most detailed.
!
!      - sll_ (not compulsary for file names)
!      - mlt_pic/lt_pic/etc
!      - 2d/4d/etc
!      - group/grid/particle/level
!      - set/get/pointer/new/delete/etc
!      - landau
!      - indices

module sll_mlt_pic_2d_particle_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_cartesian_meshes ! [[file:../meshes/sll_cartesian_meshes.F90]]

  implicit none

  !> @brief <<sll_mlt_pic_2d_particle>> One actual particle. Redefinition of \ref sll_particle_2d for multilevel and
  !> linearly-transformed particles.

  type :: sll_mlt_pic_2d_particle
     sll_int32  :: ic   ! cell index, linearly arranged
     sll_real32 :: dx
     sll_real64 :: vx
     sll_real32 :: q

     ! ALH_MCP_06_2014 Each particle knows its index in the <<remap_grid>> used in
     ! [[file:../pic_particle_initializers/sll_particle_init2D.F90]]

     sll_int32 :: k_x
     sll_int32 :: k_vx
  end type sll_mlt_pic_2d_particle

  !> @brief <<sll_mlt_pic_2d_particle_pointer>> Pointer to one [[sll_mlt_pic_2d_particle]].

  type :: sll_mlt_pic_2d_particle_pointer
     type(sll_mlt_pic_2d_particle), pointer :: p
  end type sll_mlt_pic_2d_particle_pointer

  !> @brief <<sll_mlt_pic_2d_particle_guard>> <<guard>> means a particle going out of the domain.
  !>
  !> \todo duplicates [[sll_mlt_pic_2d_particle_pointer]]?

  type :: sll_mlt_pic_2d_particle_guard
     type(sll_mlt_pic_2d_particle), pointer :: p
  end type sll_mlt_pic_2d_particle_guard

end module sll_mlt_pic_2d_particle_module
