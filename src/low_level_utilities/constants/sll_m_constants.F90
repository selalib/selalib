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

!> @ingroup constants
!> @brief
!> Fortran module where set some physical and mathematical constants.
!> @details
!> In this module, all variables must be protected
module sll_m_constants
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c, &
    sll_e_charge, &
    sll_e_mass, &
    sll_epsilon_0, &
    sll_g, &
    sll_i1, &
    sll_kb, &
    sll_kx, &
    sll_mu_0, &
    sll_pi, &
    sll_proton_mass, &
    sll_sqrt3, &
    sll_twopi

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @param PI number
sll_real64, parameter :: sll_pi = 3.1415926535897932384626433_f64 

!> @param 2*PI number
sll_real64, parameter :: sll_twopi = 2.0_f64*sll_pi

!> @param sll_kx is the fundamental mode in the x-direction.
sll_real64, parameter :: sll_kx = 2.0_f64*sll_pi

!> @param speed of light in vacuum (def) m/s  
sll_real64, parameter :: sll_c = 2.99792458e8_f64

!> @param permittivity of free space F/m      
sll_real64, parameter :: sll_epsilon_0 = 8.854187817e-12_f64

!> @param permeability of free space N/A^2      
sll_real64, parameter :: sll_mu_0 = 12.566370614e-7_f64

!> @param electron charge magnitude (49) C      
sll_real64, parameter :: sll_e_charge = 1.602176565e-19_f64

!> @param electron mass (54) kg                  
sll_real64, parameter :: sll_e_mass = 9.10938291e-31_f64

!> @param proton mass (10) kg                    
sll_real64, parameter :: sll_proton_mass = 1.672621777e-27_f64

!> @param standard grav. accel., sea level m/s^2 
sll_real64, parameter :: sll_g = 9.80665_f64

!> @param value of \f$ \sqrt{3} \f$ (needed for hexagonal meshes).
sll_real64, parameter :: sll_sqrt3 = 1.7320508075688771931766041_f64

!> @param Boltzmann constant J/K
sll_real64, parameter :: sll_kb = 1.3806488e-23_f64

!> @param Complex number i=sqrt(-1)
sll_comp64, parameter :: sll_i1 = cmplx(0.0_f64, 1.0_f64, kind=f64)

end module sll_m_constants
