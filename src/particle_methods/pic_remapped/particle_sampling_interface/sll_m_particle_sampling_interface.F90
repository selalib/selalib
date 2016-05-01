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

  use sll_m_particle_initializer, only: &
    sll_s_particle_initialize_random_landau_2d2v, &
    sll_s_particle_initialize_sobol_landau_2d2v

  use sll_m_pic_lbfr_4d_group, only: &
    sll_t_pic_lbfr_4d_group

  use sll_m_initial_density_parameters, only: &
    sll_t_initial_density_parameters, &
    sll_p_landau_density_2d2v

  implicit none

  public :: &
    sll_t_conservative_sampling_parameters, &
    sll_s_sample_particle_group, &
    sll_s_resample_particle_group, &
    sll_p_random_sampling, &
    sll_p_sobol_sampling, &
    sll_p_deterministic_sampling

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_random_sampling=0
  sll_int32, parameter :: sll_p_sobol_sampling=1
  sll_int32, parameter :: sll_p_deterministic_sampling=2

  !> type used to enforce some conservation properties in the sampling -- (there may be more than just one scalar)
  type sll_t_conservative_sampling_parameters
      sll_real64 :: total_charge      !< total charge to be conserved
    contains
      procedure :: set_total_charge
  end type sll_t_conservative_sampling_parameters


contains

  subroutine set_total_charge(self, target_total_charge)
    class(sll_t_conservative_sampling_parameters),  intent( inout )         :: self
    sll_real64,                                     intent( in ), optional  :: target_total_charge
    self%total_charge = target_total_charge
  end subroutine

  !------------------------------------------------------------------------------------------

  !> sampling interface
  subroutine sll_s_sample_particle_group( &
          particle_group,   &
          sampling_strategy, &
          initial_density_parameters, &
          conservative_sampling_parameters, &
          rank &
          )

    class(sll_c_particle_group_base),    pointer, intent( inout )         :: particle_group
    sll_int32, intent( in )                                               :: sampling_strategy
    class(sll_t_initial_density_parameters),      intent( in )            :: initial_density_parameters
    class(sll_t_conservative_sampling_parameters), intent( in ), optional :: conservative_sampling_parameters
    sll_int32, intent( in ),                                    optional  :: rank

    sll_int32, allocatable :: rnd_seed(:)
    sll_int32 :: j, ierr
    sll_int32 :: aux_rank
    sll_int32 :: rnd_seed_size
    sll_int64 :: sobol_seed

    sll_real64 :: target_total_charge
    logical    :: enforce_total_charge
    if( present(rank) )then
      aux_rank = rank
    else
      aux_rank = 0
    end if
    if( sampling_strategy == sll_p_random_sampling )then
      ! Set the seed for the random sampling
      call random_seed(size=rnd_seed_size)
      SLL_ALLOCATE(rnd_seed(rnd_seed_size), ierr)
      do j=1, rnd_seed_size
        rnd_seed(j) = (-1)**j*(100 + 15*j)*(2*aux_rank + 1)
      end do

      ! Sample position and velocity of the particles.
      ! Random sampling
      if( initial_density_parameters%f0_type == sll_p_landau_density_2d2v )then
        call sll_s_particle_initialize_random_landau_2d2v( &
            particle_group, &
            initial_density_parameters%landau_param, &
            initial_density_parameters%eta_min(1:2), &
            initial_density_parameters%domain_length(1:2), &
            initial_density_parameters%thermal_velocity, &
            rnd_seed)
      else
        SLL_ERROR("sll_s_sample_particle_group", "sampling not implemented for this initial density")
      end if

    else if( sampling_strategy == sll_p_sobol_sampling )then
      sobol_seed = int(10 + aux_rank * particle_group%n_particles, 8)

      ! Pseudo-random sampling with sobol numbers
      if( initial_density_parameters%f0_type == sll_p_landau_density_2d2v )then
        call sll_s_particle_initialize_sobol_landau_2d2v( &
            particle_group, &
            initial_density_parameters%landau_param, &
            initial_density_parameters%eta_min(1:2), &
            initial_density_parameters%domain_length(1:2), &
            initial_density_parameters%thermal_velocity, &
            sobol_seed)
      else
        SLL_ERROR("sll_s_sample_particle_group", "sampling not implemented for this initial density")
      end if

    elseif (sampling_strategy == sll_p_deterministic_sampling) then
    select type ( particle_group )

      type is ( sll_t_pic_lbfr_4d_group )
        if( present(conservative_sampling_parameters) )then
          enforce_total_charge = .true.
          target_total_charge = conservative_sampling_parameters%total_charge
        else
          enforce_total_charge = .false.    ! no charge conservation
          target_total_charge = 0._f64      ! value does not matter then
        end if
        call particle_group%resample( target_total_charge, enforce_total_charge, initial_density_parameters )

      class default
        SLL_ERROR("sll_s_sample_particle_group", "deterministic sampling interface not implemented for this type of particles")

      end select

    else
      SLL_ERROR("sll_s_sample_particle_group", "unknown strategy for the sampling interface")

    end if

  end subroutine sll_s_sample_particle_group



  !> resampling interface
  subroutine sll_s_resample_particle_group( &
          particle_group,   &
          conservative_sampling_parameters )      !< whether charge must be conserved

    class(sll_c_particle_group_base),    pointer, intent( inout )         :: particle_group
    class(sll_t_conservative_sampling_parameters), intent( in ), optional :: conservative_sampling_parameters
    sll_real64 :: target_total_charge
    logical    :: enforce_total_charge

    select type ( particle_group )

    type is ( sll_t_pic_lbfr_4d_group )
      if( present(conservative_sampling_parameters) )then
        enforce_total_charge = .true.
        target_total_charge = conservative_sampling_parameters%total_charge
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
