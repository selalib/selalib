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


module sll_particle_1d_description
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_logical_meshes

  implicit none



  type :: sll_particle_1d
     sll_real64 :: dx
     sll_real64 :: vx
     sll_real64 :: weight            !\gamma_k
     sll_real64 :: weight_const      !c_k
  end type sll_particle_1d



  type :: sll_particle_1d_guard
     type(sll_particle_1d), pointer :: p

  end type sll_particle_1d_guard


  type :: sll_particle_1d_group
    type(sll_particle_1d), allocatable, dimension(:) :: particle
    sll_real64 :: qm
  !contains
  !remove
  !add
  !
  end type sll_particle_1d_group





end module sll_particle_1d_description
