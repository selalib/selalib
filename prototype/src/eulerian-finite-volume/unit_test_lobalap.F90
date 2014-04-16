program test_lobalap
#include "selalib.h"

  use lobalap
  use map_function_module, only: set_map_function
  use sll_lobatto_poisson
  implicit none

  type(lobatto_poisson_solver)        :: solver
  type(sll_logical_mesh_2d), pointer  :: mesh
  class(sll_coordinate_transformation_2d_base), pointer :: tau

  !type(dg_field), pointer :: rho
  sll_int32 :: error
  
#define NPTS1 32
#define NPTS2 16

  ! logical mesh for space coordinates
  mesh => new_logical_mesh_2d( NPTS1, NPTS2,    & 
       eta1_min= 2.0_f64, eta1_max= 8.0_f64,         &
       eta2_min=-sll_pi, eta2_max=sll_pi)

  ! coordinate transformation associated with space coordinates
  tau => new_coordinate_transformation_2d_analytic( &
       "analytic_polar_transformation", &
       mesh, &
       polar_x1, &
       polar_x2, &
       polar_jac11, &
       polar_jac12, &
       polar_jac21, &
       polar_jac22, &
       (/ 0.0_f64 /) ) ! this particular transformation is not parametrizable

  call initialize(solver, tau)
  !rho => new_dg_field( degree, tau) 
  call solve(solver)!, rho)
  call delete(solver)


end program test_lobalap


