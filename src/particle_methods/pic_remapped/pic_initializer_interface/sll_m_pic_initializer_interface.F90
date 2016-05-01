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


module sll_m_pic_initializer_interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
! #include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_pic_lbfr_4d_group, only: &
    sll_t_pic_lbfr_4d_group

  implicit none

  public :: &
    sll_t_pic_initializer_interface

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  type sll_t_pic_initializer_interface

      !   no member fields for the moment

    contains

      procedure :: initialize_particle_group
      procedure :: set_landau_parameters
      procedure :: set_hat_f0_parameters

  end type sll_t_pic_initializer_interface


contains


  subroutine initialize_particle_group( &
          self, &
          particle_group,  &
          initial_density_identifier, &
          target_total_charge,  & !< total charge to be conserved
          enforce_total_charge  &  !< whether charge must be conserved
          )

    class(sll_t_pic_initializer_interface),     intent( inout )        :: self
    class(sll_c_particle_group_base), pointer,  intent( inout )        :: particle_group
    sll_int32,                                  intent( in    )        :: initial_density_identifier
    sll_real64,                                 intent( in ), optional :: target_total_charge
    logical,                                    intent( in ), optional :: enforce_total_charge

    sll_real64 :: aux_target_total_charge
    logical    :: aux_enforce_total_charge

    select type ( particle_group )

    type is ( sll_t_pic_lbfr_4d_group )
      if( present(target_total_charge) )then
        SLL_ASSERT( present(enforce_total_charge) )
        aux_enforce_total_charge = enforce_total_charge
        aux_target_total_charge = target_total_charge
      else
        aux_enforce_total_charge = .false.      ! default: does not enforce charge conservation
        aux_target_total_charge = 0._f64        ! value does not matter then
      end if
      call particle_group%initializer( &
          initial_density_identifier,   &
          aux_target_total_charge,  &
          aux_enforce_total_charge  &
      )

    class default
      SLL_ERROR("initialize_particle_group", "initializer procedure not implemented for this type of particle group")

    end select

  end subroutine initialize_particle_group

  subroutine set_landau_parameters( &
          self, &
          particle_group,  &
          thermal_speed, alpha, k_landau &
          )

    class(sll_t_pic_initializer_interface),     intent( inout )        :: self
    class(sll_c_particle_group_base), pointer,  intent( inout )        :: particle_group
    sll_real64,                                 intent( in    )        :: thermal_speed
    sll_real64,                                 intent( in    )        :: alpha
    sll_real64,                                 intent( in    )        :: k_landau

    select type ( particle_group )

    type is ( sll_t_pic_lbfr_4d_group )
      call particle_group%set_landau_parameters( thermal_speed, alpha, k_landau )

    class default
      SLL_ERROR("set_landau_parameters", "procedure not implemented for this type of particle group")

    end select

  end subroutine set_landau_parameters

  subroutine set_hat_f0_parameters( &
          self, &
          particle_group,  &
          x0, y0, vx0, vy0, r_x, r_y, r_vx, r_vy, basis_height, shift &
          )

    class(sll_t_pic_initializer_interface),       intent( inout )        :: self
    class(sll_c_particle_group_base), pointer,    intent( inout )        :: particle_group
    sll_real64,                                   intent( in    )        :: x0
    sll_real64,                                   intent( in    )        :: y0
    sll_real64,                                   intent( in    )        :: vx0
    sll_real64,                                   intent( in    )        :: vy0
    sll_real64,                                   intent( in    )        :: r_x
    sll_real64,                                   intent( in    )        :: r_y
    sll_real64,                                   intent( in    )        :: r_vx
    sll_real64,                                   intent( in    )        :: r_vy
    sll_real64,                                   intent( in    )        :: basis_height
    sll_real64,                                   intent( in    )        :: shift

    select type ( particle_group )

    type is ( sll_t_pic_lbfr_4d_group )
      call particle_group%set_hat_f0_parameters( x0, y0, vx0, vy0, r_x, r_y, r_vx, r_vy, basis_height, shift )

    class default
      SLL_ERROR("set_hat_f0_parameters", "procedure not implemented for this type of particle group")

    end select

  end subroutine set_hat_f0_parameters

end module  sll_m_pic_initializer_interface
