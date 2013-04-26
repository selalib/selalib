#ifndef _poisson_solvers_h_
#define _poisson_solvers_h_

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

#define GET_POISSON_MESH_DESCRIPTOR( p )     p%descriptor

#define GET_POISSON_RHS( p )       p%rhs
#define GET_POISSON_SOLUTION( p )  p%sol
 
#define SET_POISSON_MESH_DESCRIPTOR( p, m )     p%descriptor => m
#define SET_POISSON_RHS( p, f )    p%rhs => f

use sll_poisson_1d_periodic
use sll_poisson_2d_periodic
use sll_poisson_2d_polar

#endif
