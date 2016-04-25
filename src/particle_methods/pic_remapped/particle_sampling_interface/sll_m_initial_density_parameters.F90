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
!> @brief Type for initial densities

module sll_m_initial_density_parameters

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_t_initial_density_parameters, &
    sll_p_landau_density_2d2v, &
    sll_p_hat_density_2d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_landau_density_2d2v=780
  sll_int32, parameter :: sll_p_hat_density_2d2v=781      !< hat function density for testing purposes

  !> type describing the initial density to be sampled
  type sll_t_initial_density_parameters
      sll_int32     :: f0_type

      sll_real64    :: eta_min(3)             !< lower bound of the domain with 3 possible dimensions
      sll_real64    :: domain_length(3)       !< length of the domain

      !> parameters for the landau density
      sll_real64    :: landau_param(2)        !< parameter defining the perturbation: landau_param(1)*cos(landau_param(2)*x1)
      sll_real64    :: thermal_velocity(3)    !< Value of the thermal velocity along each dimension.

      !> parameters for the hat density (for testing purposes)
      sll_real64    :: hat_f0_x0
      sll_real64    :: hat_f0_y0
      sll_real64    :: hat_f0_vx0
      sll_real64    :: hat_f0_vy0
      sll_real64    :: hat_f0_r_x
      sll_real64    :: hat_f0_r_y
      sll_real64    :: hat_f0_r_vx
      sll_real64    :: hat_f0_r_vy
      sll_real64    :: hat_f0_basis_height
      sll_real64    :: hat_f0_shift

    contains

      procedure :: set_landau_2d2v_parameters
      procedure :: set_hat_density_2d2v_parameters

  end type sll_t_initial_density_parameters

contains


  !> define the initial density for a landau 2d2v test case
  subroutine set_landau_2d2v_parameters(  &
        self, &
        landau_param, &
        thermal_velocity,  &
        eta1_min, &
        eta1_max, &
        eta2_min, &
        eta2_max )
    class(sll_t_initial_density_parameters),  intent( inout )         :: self
    sll_real64,   intent( in ) :: landau_param(2)
    sll_real64,   intent( in ) :: thermal_velocity(2)
    sll_real64,   intent( in ) :: eta1_min
    sll_real64,   intent( in ) :: eta1_max
    sll_real64,   intent( in ) :: eta2_min
    sll_real64,   intent( in ) :: eta2_max

    self%f0_type = sll_p_landau_density_2d2v
    self%landau_param = landau_param
    self%thermal_velocity(1:2) = thermal_velocity
    self%eta_min(1) = eta1_min
    self%eta_min(2) = eta2_min
    self%eta_min(3) = 0._f64

    self%domain_length(1) = eta1_max - eta1_min
    self%domain_length(2) = eta2_max - eta2_min
    self%domain_length(3) = 0._f64

  end subroutine set_landau_2d2v_parameters


  !> define the initial hat density for testing
  subroutine set_hat_density_2d2v_parameters(   &
        self, &
        x0, &
        y0, &
        vx0, &
        vy0, &
        r_x, &
        r_y, &
        r_vx, &
        r_vy, &
        basis_height, &
        shift )
    class(sll_t_initial_density_parameters), intent(inout) :: self
    sll_real64,                     intent(in)      :: x0
    sll_real64,                     intent(in)      :: y0
    sll_real64,                     intent(in)      :: vx0
    sll_real64,                     intent(in)      :: vy0
    sll_real64,                     intent(in)      :: r_x
    sll_real64,                     intent(in)      :: r_y
    sll_real64,                     intent(in)      :: r_vx
    sll_real64,                     intent(in)      :: r_vy
    sll_real64,                     intent(in)      :: basis_height
    sll_real64,                     intent(in)      :: shift

    self%hat_f0_x0 = x0
    self%hat_f0_y0 = y0
    self%hat_f0_vx0 = vx0
    self%hat_f0_vy0 = vy0
    self%hat_f0_r_x = r_x
    self%hat_f0_r_y = r_y
    self%hat_f0_r_vx = r_vx
    self%hat_f0_r_vy = r_vy
    self%hat_f0_basis_height = basis_height
    self%hat_f0_shift = shift

  end subroutine set_hat_density_2d2v_parameters


end module  sll_m_initial_density_parameters
