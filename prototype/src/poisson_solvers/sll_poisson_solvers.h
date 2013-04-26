#ifndef _poisson_solvers_h_
#define _poisson_solvers_h_

#define GET_POISSON_MESH_DESCRIPTOR( p )     p%descriptor

#define GET_POISSON_RHS( p )       p%rhs
#define GET_POISSON_SOLUTION( p )  p%sol
 
#define SET_POISSON_MESH_DESCRIPTOR( p, m )     p%descriptor => m
#define SET_POISSON_RHS( p, f )    p%rhs => f

use sll_poisson_1d_periodic
use sll_poisson_2d_periodic
use sll_poisson_2d_polar

#endif
