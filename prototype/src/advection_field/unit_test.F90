program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
  use sll_constants
  use sll_advection_field
!  use sll_module_mapped_meshes_2d
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_module_mapped_meshes_1d
  use sll_cubic_spline_interpolator_1d
  use sll_common_coordinate_transformations
  implicit none
 
  sll_int32 :: nc_eta1, nc_eta2
  !type(sll_mapped_mesh_2d_analytic),target     :: mesh2d
  !class(sll_mapped_mesh_2d_base), pointer      :: pm2d
  type(sll_logical_mesh_2d), pointer :: mesh2d
  class(sll_coordinate_transformation_2d_base), pointer   :: pm2d
  class(sll_mapped_mesh_1d_base), pointer      :: pm1d
  !class(scalar_field_2d_initializer_base), pointer    :: pfinit
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
  type(hamiltonian_advection_field_2d)         :: adv_field
  type(sll_mapped_mesh_1d_analytic), target    :: mesh1d
  type(scalar_field_1d)                        :: phi_self
  character(len=32)                            :: name = 'adv_field'

  type(cubic_spline_1d_interpolator), target   :: interp_eta1
  type(cubic_spline_1d_interpolator), target   :: interp_eta2
  sll_real64, dimension(4)                     :: affine_map_params
  sll_real64, dimension(2)                     :: linear_map_params
  sll_int32                                    :: ix
  sll_int32, parameter                         :: NODE_CENTERED_FIELD_COPY = 0

  nc_eta1 = 100
  nc_eta2 = 100

  affine_map_params(:) = (/-1.0_f64, 1.0_f64, -1.0_f64, 1.0_f64/)

  print*, 'initialization of mesh'
  mesh2d => new_logical_mesh_2d( &
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

  print*, 'define 1d mesh in x for potential'
  call mesh1d%initialize( &
       "mesh1d",  &
       nc_eta1+1, &
       linear_map_f,      &
       linear_map_jac_f,  &
       linear_map_params )
  pm1d => mesh1d

  print*, 'initialization of 1D potential'
  call initialize_scalar_field_1d( &
       phi_self, &
       "phi_self", &
       pm1d, &
       NODE_CENTERED_FIELD_COPY)

  do ix = 1, nc_eta1
     phi_self%data(ix) = pm1d%x1_at_node(ix)**2
  end do

  call compute_hamiltonian(adv_field, phi_self)

  call write_scalar_field_2d(adv_field) 

    
  print *, 'Successful, exiting program.'
  
end program unit_test
