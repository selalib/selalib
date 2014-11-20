program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
  use sll_constants
  use sll_advection_field
  use sll_cartesian_meshes
  use sll_module_coordinate_transformations_2d
  use sll_module_cubic_spline_interpolator_1d
  use sll_common_coordinate_transformations
  implicit none
 
  sll_int32 :: nc_eta1, nc_eta2
  type(sll_cartesian_mesh_2d), pointer :: mesh2d
  class(sll_coordinate_transformation_2d_base), pointer   :: pm2d
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
  type(hamiltonian_advection_field_2d)         :: adv_field
  character(len=32)                            :: name = 'adv_field'

  type(sll_cubic_spline_interpolator_1d), target   :: interp_eta1
  type(sll_cubic_spline_interpolator_1d), target   :: interp_eta2
  sll_real64, dimension(5)                     :: affine_map_params
  sll_real64, dimension(2)                     :: linear_map_params
  sll_int32                                    :: ix
  sll_int32, parameter                         :: NODE_CENTERED_FIELD_COPY = 0

  nc_eta1 = 100
  nc_eta2 = 100

  affine_map_params(:) = (/-1.0_f64, 1.0_f64, -1.0_f64, 1.0_f64, 0.0_f64/)

  print*, 'initialization of mesh'
  mesh2d => new_cartesian_mesh_2d( &
       nc_eta1, &
       nc_eta2  &
   )
  pm2d => new_coordinate_transformation_2d_analytic( &
       "mesh2d_cart",      &
       mesh2d,             &
       affine_x1, &
       affine_x2, &
       affine_jac11, &
       affine_jac12, &
       affine_jac21, &
       affine_jac22, &
       affine_map_params )

  ! Set up the interpolators for the field
  call interp_eta1%initialize( nc_eta1+1, 0.0_f64, 1.0_f64, SLL_PERIODIC )
  call interp_eta2%initialize( nc_eta2+1, 0.0_f64, 1.0_f64, SLL_PERIODIC )
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2

  print*, 'initialization of advection field'

  call initialize_advection_field_2d( &
       adv_field, &
       1.0_f64, &
       1.0_f64, &
       name, &
       pm2d, &
       NODE_CENTERED_FIELD_COPY, &
       eta1_interpolator=interp_eta1_ptr, &
       eta2_interpolator=interp_eta2_ptr )
  
  linear_map_params(:) = (/-1.0_f64,1.0_f64/)


  call write_scalar_field_2d(adv_field) 

  print *, 'Successful, exiting program.'
  
end program unit_test
