#ifndef _SELALIB_MPI_H_
#define _SELALIB_MPI_H_

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


#include "selalib.h"

use sll_collective
use sll_remapper
use sll_gnuplot_parallel
use sll_poisson_2d_periodic_cartesian_par

#define MPI_MASTER 0

#endif

