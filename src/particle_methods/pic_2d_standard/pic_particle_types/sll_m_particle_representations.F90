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


module sll_m_particle_representations
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_particle_2d, &
    sll_particle_2d_guard, &
    sll_particle_2d_guard_ptr, &
    sll_particle_4d, &
    sll_particle_4d_guard, &
    sll_particle_4d_guard_ptr

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type :: sll_particle_4d
     sll_int32  :: ic   ! cell index, linearly arranged
     sll_real32 :: dx!  sll_real64 :: dx!
     sll_real32 :: dy!  sll_real64 :: dy
     sll_real64 :: vx
     sll_real64 :: vy
     sll_real32 :: q!   sll_real64 :: q!      
  end type sll_particle_4d

  type :: sll_particle_4d_guard
     type(sll_particle_4d), pointer :: p
  end type sll_particle_4d_guard

  type :: sll_particle_4d_guard_ptr
     type(sll_particle_4d_guard), dimension(:), pointer :: g_list
  end type sll_particle_4d_guard_ptr


! ------------------------------
!  for the GUIDING CENTER model
! ------------------------------
  type :: sll_particle_2d
     sll_int32  :: ic   ! cell index, linearly arranged
     sll_real32 :: dx
     sll_real32 :: dy
     sll_real32 :: q
  end type sll_particle_2d

  type :: sll_particle_2d_guard
     type(sll_particle_2d), pointer :: p
  end type sll_particle_2d_guard

  type :: sll_particle_2d_guard_ptr
     type(sll_particle_2d_guard), dimension(:), pointer :: g_list
  end type sll_particle_2d_guard_ptr


!contains

end module sll_m_particle_representations
