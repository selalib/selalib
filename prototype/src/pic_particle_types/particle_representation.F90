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


module sll_particle_representations
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  type :: sll_particle_2d
     sll_int32 :: ic   ! cell index, linearly arranged
     sll_real32 :: dx
     sll_real32 :: dy
     sll_real64 :: vx
     sll_real64 :: vy
     sll_real32 :: q
  end type sll_particle_2d

  type :: sll_particle_2d_guard
     type(sll_particle_2d), pointer :: p
  end type sll_particle_2d_guard

contains

  subroutine initialize_particle_2d( &
       p,  &
       ic, &
       dx, &
       dy, &
       vx, &
       vy, &
       q )

    type(sll_particle_2d), intent(out) :: p
    sll_int32,  intent(in) :: ic   ! cell index, linearly arranged
    sll_real32, intent(in) :: dx
    sll_real32, intent(in) :: dy
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real32, intent(in) :: q

    p%ic = ic
    p%dx = dx
    p%dy = dy
    p%vx = vx
    p%vy = vy
    p%q  = q
  end subroutine initialize_particle_2d

end module sll_particle_representations
