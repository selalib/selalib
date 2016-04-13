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


module sll_m_pic_resamplers

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
! #include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_remapped_pic_base, only: &
    sll_c_remapped_particle_group

  use sll_m_pic_lbfr_4d_group, only: &
    sll_t_pic_lbfr_4d_group

  implicit none

  public :: &
    sll_t_pic_4d_resampler

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  type sll_t_pic_4d_resampler

      !   no member fields for the moment

    contains

      procedure :: resample_particle_group

  end type sll_t_pic_4d_resampler


contains


  subroutine resample_particle_group( &
          self, &
          particle_group,   &
          target_total_charge,  &     !< total charge to be conserved
          enforce_total_charge )      !< whether charge must be conserved

    class(sll_t_pic_4d_resampler),                  intent( inout )         :: self
    class(sll_c_remapped_particle_group),  pointer, intent( inout )         :: particle_group
    sll_real64,                                     intent( in ), optional  :: target_total_charge
    logical,                                        intent( in ), optional  :: enforce_total_charge
    sll_real64  :: aux_target_total_charge
    logical     :: aux_enforce_total_charge

    select type ( particle_group )
    type is ( sll_t_pic_lbfr_4d_group )
      if( present(target_total_charge) )then
        SLL_ASSERT( present(enforce_total_charge) )
        aux_enforce_total_charge = enforce_total_charge
        aux_target_total_charge = target_total_charge
      else
        aux_enforce_total_charge = .false.    ! default: does not enforce charge conservation
        aux_target_total_charge = 0._f64           ! value does not matter then
      end if
      call particle_group%remap( aux_target_total_charge, aux_enforce_total_charge )
    class default
      SLL_ERROR("pic_4d_resampling", "resampling procedure should not be called for this type of particle group")
    end select

  end subroutine resample_particle_group


end module  sll_m_pic_resamplers
