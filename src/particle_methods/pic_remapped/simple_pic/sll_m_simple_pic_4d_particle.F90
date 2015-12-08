!> @file simple_pic_4d_particle.F90

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

!> @details Mainly useful to compare against more advanced models. A group of sll_simple_pic_4d_particle
!> is a @ref sll_simple_pic_4d_group

! [[sll_simple_pic_4d_particle]]
! [[selalib:src/particle_methods/particle_types/simple_pic_4d_group.F90::sll_simple_pic_4d_group]]

module sll_m_simple_pic_4d_particle

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_simple_pic_4d_particle

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! <<sll_simple_pic_4d_particle>>

  !> @ingroup particle_methods

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
     sll_real64 :: vx

     !> speed along Y axis
     sll_real64 :: vy

     sll_real32 :: weight

  end type sll_simple_pic_4d_particle

end module sll_m_simple_pic_4d_particle
