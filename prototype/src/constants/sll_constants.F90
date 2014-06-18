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

!> @brief
!> Fortran module where set some physical and mathematical constants.
!> @details
!> This file should contain all mathematical and physical
!> constants to be used library-wide.
!>
!> <b> How to use it </b>
!> @code
!> #include "sll_constants.h"
!> @endcode
module sll_constants
#include "sll_working_precision.h"

implicit none

!> @param PI number
sll_real64, parameter :: sll_pi = 3.1415926535897932384626433_f64

!> @param sll_kx is the fundamental mode in the x-direction.
sll_real64            :: sll_kx = 2*sll_pi

!> @param speed of light in vacuum (def) m/s
sll_real64, parameter :: sll_c = 2.99792458D8

!> @param permittivity of free space F/m
sll_real64, parameter :: sll_epsilon_0 = 8.854187817D-12

!> @param permeability of free space N/A^2
sll_real64, parameter :: sll_mu_0 = 12.566370614D-7

!> @param electron charge magnitude (49) C
sll_real64, parameter :: sll_e_charge = 1.602176565D-19

!> @param electron mass (54) kg
sll_real64, parameter :: sll_e_mass = 9.10938291D-31

!> @param proton mass (10) kg
sll_real64, parameter :: sll_proton_mass = 1.672621777D-27

!> @param standard grav. accel., sea level m/s^2
sll_real64, parameter :: sll_g = 9.80665D0

!> @param Boltzmann constant J/K
sll_real64, parameter :: sll_kb = 1.3806488D-23

!> @param Complex number i=sqrt(-1)
sll_comp64, parameter :: sll_i1 = (0D0, 1D0)

end module sll_constants
