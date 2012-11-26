program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
  use numeric_constants
  use sll_advection_field
  use sll_module_mapped_meshes_2d
  use sll_module_mapped_meshes_1d
  use sll_cubic_spline_interpolator_1d
  use geometry_functions
  implicit none
 
  sll_int32 :: nc_eta1, nc_eta2
  type(sll_mapped_mesh_2d_analytic),target     :: mesh2d
#ifdef STDF95
  type(sll_mapped_mesh_2d_analytic), pointer      :: pm2d
  type(sll_mapped_mesh_1d_analytic), pointer      :: pm1d
  type(scalar_field_2d_initializer_base), pointer    :: pfinit
  type(cubic_spline_1d_interpolator), pointer :: interp_eta1_ptr
  type(cubic_spline_1d_interpolator), pointer :: interp_eta2_ptr
#else
  class(sll_mapped_mesh_2d_base), pointer      :: pm2d
  class(sll_mapped_mesh_1d_base), pointer      :: pm1d
  class(scalar_field_2d_initializer_base), pointer    :: pfinit
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
#endif
  type(hamiltonian_advection_field_2d)         :: adv_field
  type(sll_mapped_mesh_1d_analytic), target    :: mesh1d
  type(scalar_field_1d)                        :: phi_self
  character(len=32)                            :: name = 'adv_field'
  character(len=4)                             :: cstep

  type(cubic_spline_1d_interpolator), target  :: interp_eta1
  type(cubic_spline_1d_interpolator), target  :: interp_eta2


  sll_int32                                    :: ierr, istep
  sll_int32                                    :: ix, iv, nnode_x1, nnode_v1
  sll_int32, parameter                         :: NODE_CENTERED_FIELD_COPY = 0

  nc_eta1 = 100
  nc_eta2 = 100

  print*, 'initialization of mesh'

#ifdef STDF95
  call initialize_mesh_2d_analytic( &
       mesh2d, &
#else
  call mesh2d%initialize( &
#endif
       "mesh2d",  &
       nc_eta1+1, &
       nc_eta2+1, &
       affine_x1, &
       affine_x2, &
       affine_jac11, &
       affine_jac12, &
       affine_jac21, &
       affine_jac22 )
  pm2d => mesh2d

  ! Set up the interpolators for the field
#ifdef STDF95
  call cubic_spline_1d_interpolator_initialize(interp_eta1, nc_eta1+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  call cubic_spline_1d_interpolator_initialize(interp_eta2, nc_eta2+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
#else
  call interp_eta1%initialize( nc_eta1+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  call interp_eta2%initialize( nc_eta2+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
#endif
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
  
  print*, 'define 1d mesh in x for potential'
#ifdef STDF95
  call initialize_mesh_1d_analytic( &
       mesh1d, &
#else
  call mesh1d%initialize( &
#endif
       "mesh1d",  &
       nc_eta1+1, &
       linear_map_f,        &
       linear_map_jac_f )
  pm1d => mesh1d

  print*, 'initialization of 1D potential'
  call initialize_scalar_field_1d( &
       phi_self, &
       "phi_self", &
       pm1d, &
       NODE_CENTERED_FIELD_COPY)

  do ix = 1, nc_eta1
#ifdef STDF95
     phi_self%data(ix) = x1_node_analytic_1d(pm1d, ix)**2
#else
     phi_self%data(ix) = pm1d%x1_at_node(ix)**2
#endif
  end do

  call compute_hamiltonian(adv_field, phi_self)
#ifdef STDF95
  call write_scalar_field_2d(adv_field%extend_type) 
#else
  call write_scalar_field_2d(adv_field) 
#endif

    
  print *, 'Successful, exiting program.'
  
end program unit_test
