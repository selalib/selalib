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
      sll_p_c, &
      sll_p_charge, &
      sll_p_mass, &
      sll_p_epsilon_0, &
      sll_p_fourpi, &
      sll_p_g, &
      sll_p_i0, &
      sll_p_i1, &
      sll_p_kb, &
      sll_p_kx, &
      sll_p_mu_0, &
      sll_p_pi, &
      sll_p_proton_mass, &
      sll_p_sqrt3, &
      sll_p_twopi

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @param PI number
!sll_real64, parameter :: sll_p_pi = 3.1415926535897932384626433_f64
   sll_real64, parameter :: sll_p_pi = 3.141592653589793_f64

!> @param 2*PI number
   sll_real64, parameter :: sll_p_twopi = 2.0_f64*sll_p_pi

!> @param 4*PI number
   sll_real64, parameter :: sll_p_fourpi = 4.0_f64*sll_p_pi

!> @param sll_p_kx is the fundamental mode in the x-direction.
   sll_real64, parameter :: sll_p_kx = 2.0_f64*sll_p_pi

!> @param speed of light in vacuum (def) m/s
   sll_real64, parameter :: sll_p_c = 2.99792458e8_f64

!> @param permittivity of free space F/m
   sll_real64, parameter :: sll_p_epsilon_0 = 8.854187817e-12_f64

!> @param permeability of free space N/A^2
   sll_real64, parameter :: sll_p_mu_0 = 12.566370614e-7_f64

!> @param electron charge magnitude (49) C
   sll_real64, parameter :: sll_p_charge = 1.602176565e-19_f64

!> @param electron mass (54) kg
   sll_real64, parameter :: sll_p_mass = 9.10938291e-31_f64

!> @param proton mass (10) kg
   sll_real64, parameter :: sll_p_proton_mass = 1.672621777e-27_f64

!> @param standard grav. accel., sea level m/s^2
   sll_real64, parameter :: sll_p_g = 9.80665_f64

!> @param value of \f$ \sqrt{3} \f$ (needed for hexagonal meshes).
!sll_real64, parameter :: sll_p_sqrt3 = 1.7320508075688771931766041_f64
   sll_real64, parameter :: sll_p_sqrt3 = 1.7320508075688772_f64

!> @param Boltzmann constant J/K
   sll_real64, parameter :: sll_p_kb = 1.3806488e-23_f64

!> @param Complex number i=sqrt(-1)
   sll_comp64, parameter :: sll_p_i1 = cmplx(0.0_f64, 1.0_f64, kind=f64)

!> @param Complex number i=sqrt(-1)
   sll_comp64, parameter :: sll_p_i0 = (0.0_f64, 0.0_f64)

end module sll_m_constants
