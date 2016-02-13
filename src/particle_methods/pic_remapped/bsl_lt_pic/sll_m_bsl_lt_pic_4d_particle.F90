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

!> @brief Particle model for bsl_lt_pic method

module sll_m_bsl_lt_pic_4d_particle

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_t_bsl_lt_pic_4d_particle

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! <<sll_t_bsl_lt_pic_4d_particle>>
  type :: sll_t_bsl_lt_pic_4d_particle
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

    sll_real64 :: weight

    !> @name Particle indices of logical neighbors (= the neighbors on the cartesian remapping grid)
    !> @{
    sll_int32  :: ngb_xleft_index
    sll_int32  :: ngb_xright_index
    sll_int32  :: ngb_yleft_index
    sll_int32  :: ngb_yright_index
    sll_int32  :: ngb_vxleft_index
    sll_int32  :: ngb_vxright_index
    sll_int32  :: ngb_vyleft_index
    sll_int32  :: ngb_vyright_index
    !> @}

  end type sll_t_bsl_lt_pic_4d_particle

end module sll_m_bsl_lt_pic_4d_particle
