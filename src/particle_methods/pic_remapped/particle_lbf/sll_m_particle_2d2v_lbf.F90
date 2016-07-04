!**************************************************************
!  Copyright INRIA
!  Authors : 
!     MCP ALH
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

!> @brief particles used for linearized-backward-flow (lbf) resamplings on cartesian grids

module sll_m_particle_2d2v_lbf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_t_particle_2d2v_lbf

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! <<sll_t_particle_2d2v_lbf>>
  type :: sll_t_particle_2d2v_lbf

    sll_real64, dimension(4) :: eta
    sll_real64, dimension(1) :: weights   !> charge = species_charge * weights(1) -- only one weight for the moment

    !> @name indices of logical neighbors (= the neighbors on the cartesian grid)
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

  end type sll_t_particle_2d2v_lbf

end module sll_m_particle_2d2v_lbf
