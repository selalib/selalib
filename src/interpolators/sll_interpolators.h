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

use sll_boundary_condition_descriptors
use sll_module_interpolators_1d_base
use sll_module_cubic_spline_interpolator_1d
use sll_module_cubic_spline_interpolator_1d_nonuniform
use sll_module_arbitrary_degree_spline_interpolator_1d
use sll_module_interpolators_2d_base
use sll_module_cubic_spline_interpolator_2d
use sll_module_arbitrary_degree_spline_interpolator_2d
!use sll_module_bspline_interpolator_1d
!use sll_module_bspline_interpolator_2d
use sll_module_periodic_interpolator_1d
use sll_module_lagrange_interpolator_1d
use periodic_interp_module
use sll_lagrange_interpolation

#include "sll_splines.h"

#endif
 
