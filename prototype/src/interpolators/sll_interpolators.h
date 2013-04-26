#ifndef _sll_interpolators_h_
#define _sll_interpolators_h_

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

#ifndef STDF95
use sll_module_interpolators_1d_base
#endif
use sll_cubic_spline_interpolator_1d
use sll_quintic_spline_interpolator_1d
use sll_odd_degree_spline_interpolator_1d

#ifndef STDF95
use sll_module_interpolators_2d_base
#endif
use sll_cubic_spline_interpolator_2d

#include "sll_splines.h"

#endif
 
