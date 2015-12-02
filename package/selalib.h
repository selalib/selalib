#ifndef _SELALIB_H_
#define _SELALIB_H_

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

#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_m_cartesian_meshes
use sll_m_working_precision
use sll_m_utilities
use sll_m_constants
use sll_m_gnuplot
use sll_m_cubic_splines
use sll_m_cubic_spline_interpolator_1d
use sll_m_cubic_spline_interpolator_2d
use sll_m_interpolators_1d_base
use sll_m_cubic_splines
use sll_m_coordinate_transformation_2d_base
use sll_m_coordinate_transformations_2d
use sll_m_common_coordinate_transformations
use sll_m_poisson_1d_periodic
use sll_m_poisson_2d_periodic_fftpack

#define poisson_2d_periodic poisson_2d_periodic_fftpack


#endif
