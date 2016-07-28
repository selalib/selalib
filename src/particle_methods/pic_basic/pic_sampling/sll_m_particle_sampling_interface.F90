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

!> @ingroup particle_methods
!> @author MCP
!> @brief Interface routines for sampling and resampling particle groups.

module sll_m_particle_sampling_interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_sampling, only: &
       sll_t_particle_sampling

  use sll_m_particle_group_2d2v_lbf, only: &
    sll_t_particle_group_2d2v_lbf

  use sll_m_initial_distribution, only : &
       sll_c_distribution_params

  implicit none

  public :: &
    sll_t_conservative_sampling_params, &
    sll_s_sample_particle_group, &
    sll_s_resample_particle_group

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> type used to enforce some conservation properties in the sampling -- (there may be more than just one scalar)
  type sll_t_conservative_sampling_params
      sll_real64 :: total_charge      !< total charge to be conserved
    contains
      procedure :: set_total_charge
  end type sll_t_conservative_sampling_params


contains

  subroutine set_total_charge(self, target_total_charge)
    class(sll_t_conservative_sampling_params),  intent( inout )         :: self
    sll_real64,                                     intent( in ), optional  :: target_total_charge
    self%total_charge = target_total_charge
  end subroutine

  !------------------------------------------------------------------------------------------

  !> sampling interface
  subroutine sll_s_sample_particle_group( &
          particle_group,   &
          initial_distribution_params, &
          random_sampler, &
          conservative_sampling_params, &
          xmin, &
          Lx )

    class(sll_c_particle_group_base), pointer,  intent( inout )        :: particle_group
    class(sll_c_distribution_params),           intent( in )           :: initial_distribution_params
    type(sll_t_particle_sampling),              intent( inout ), optional :: random_sampler   !< must have been initialized before
    class(sll_t_conservative_sampling_params),  intent( in ), optional :: conservative_sampling_params
    sll_real64,                                 intent( in ), optional :: xmin(:)  !< lower bound of the domain
    sll_real64,                                 intent( in ), optional :: Lx(:)    !< length of the domain.

    sll_real64 :: target_total_charge
    logical    :: enforce_total_charge

    select type ( particle_group )

      type is ( sll_t_particle_group_2d2v_lbf )

        !> default sampling strategy for lbf particles is type-bound deterministic, no use of (optional) random_sampler object
        if( present(conservative_sampling_params) )then
          enforce_total_charge = .true.
          target_total_charge = conservative_sampling_params%total_charge
        else
          enforce_total_charge = .false.    ! no charge conservation
          target_total_charge = 0._f64      ! value does not matter then
        end if
        call particle_group%sample( target_total_charge, enforce_total_charge, initial_distribution_params )

      class default

        !> default sample using the (optional) random_sampler object
        SLL_ASSERT( present( random_sampler ) )
        SLL_ASSERT( present( xmin ) )
        SLL_ASSERT( present( Lx ) )
        call random_sampler%sample( particle_group, initial_distribution_params, xmin, Lx )

      end select

  end subroutine sll_s_sample_particle_group



  !> resampling interface
  subroutine sll_s_resample_particle_group( &
          particle_group,   &
          conservative_sampling_params )      !< whether charge must be conserved

    class(sll_c_particle_group_base),    pointer, intent( inout )         :: particle_group
    class(sll_t_conservative_sampling_params), intent( in ), optional :: conservative_sampling_params
    sll_real64 :: target_total_charge
    logical    :: enforce_total_charge

    select type ( particle_group )

    type is ( sll_t_particle_group_2d2v_lbf )
      if( present(conservative_sampling_params) )then
        enforce_total_charge = .true.
        target_total_charge = conservative_sampling_params%total_charge
      else
        enforce_total_charge = .false.    ! no charge conservation
        target_total_charge = 0._f64      ! value does not matter then
      end if
      call particle_group%resample( target_total_charge, enforce_total_charge )

    class default
      SLL_ERROR("sll_s_resample_particle_group", "resampling interface not implemented for this type of particle group")

    end select

  end subroutine sll_s_resample_particle_group



end module  sll_m_particle_sampling_interface
