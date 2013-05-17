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

module sll_constants
#include "sll_working_precision.h"

  implicit none


  sll_real64, parameter :: sll_pi = 3.1415926535897932384626433_f64
  ! sll_kx is the fundamental mode in the x-direction. 
  ! It should be set at mesh initialization
  sll_real64            :: sll_kx = 2*sll_pi 


  sll_real64, parameter :: sll_c = 2.99792458D8            ! speed of light in vacuum (def) m/s  
  sll_real64, parameter :: sll_epsilon_0 = 8.854187817D-12 ! permittivity of free space F/m      
  sll_real64, parameter :: sll_mu_0 = 12.566370614D-7      ! permeability of free space N/A^2      
  sll_real64, parameter :: sll_e_charge = 1.60217733D-19   ! electron charge magnitude (49) C      
  sll_real64, parameter :: sll_e_mass = 9.1093897D-31      ! electron mass (54) kg                  
  sll_real64, parameter :: sll_proton_mass = 1.6726231D-27 ! proton mass (10) kg                    
  sll_real64, parameter :: sll_g = 9.80665D0               ! standard grav. accel., sea level m/s^2 

end module sll_constants
