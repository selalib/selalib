program unit_test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_common_coordinate_transformations, only: &
    sll_f_sinprod_jac11, &
    sll_f_sinprod_jac12, &
    sll_f_sinprod_jac21, &
    sll_f_sinprod_jac22, &
    sll_f_sinprod_x1, &
    sll_f_sinprod_x2

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_t_cubic_spline_interpolator_1d

  use sll_m_distribution_function, only: &
    sll_s_initialize_distribution_function_2d, &
    sll_t_distribution_function_2d

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_landau_2d_initializer, only: &
    sll_t_init_landau_2d

  use sll_m_scalar_field_2d_old, only: &
    sll_s_write_scalar_field_2d

  use sll_m_scalar_field_initializers_base, only: &
    sll_p_cell_centered_field, &
    sll_c_scalar_field_2d_initializer_base

  use sll_m_utilities, only: &
    sll_s_int2string

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  sll_int32 :: nc_eta1, nc_eta2
  class(sll_c_coordinate_transformation_2d_base), pointer   :: m
  type(sll_t_cartesian_mesh_2d), pointer :: m_log
  class(sll_c_scalar_field_2d_initializer_base), pointer    :: p_init_f
  class(sll_c_interpolator_1d), pointer :: interp_eta1_ptr
  class(sll_c_interpolator_1d), pointer :: interp_eta2_ptr
  type(sll_t_distribution_function_2d)   :: df 
  character(32)  :: name = 'dist_func'
  character(len=4) :: cstep
  type(sll_t_init_landau_2d), target :: init_landau
  type(sll_t_cubic_spline_interpolator_1d), target  :: interp_eta1
  type(sll_t_cubic_spline_interpolator_1d), target  :: interp_eta2

  sll_int32 :: istep

  nc_eta1 = 100
  nc_eta2 = 100

  print*, 'initialization of mesh'
  
 m_log => sll_f_new_cartesian_mesh_2d( &
       nc_eta1, &
       nc_eta2  &
   )
  m => sll_f_new_coordinate_transformation_2d_analytic( &
       "mesh2d_coll",      &
       m_log,             &
       sll_f_sinprod_x1, &
       sll_f_sinprod_x2, &
       sll_f_sinprod_jac11, &
       sll_f_sinprod_jac12, &
       sll_f_sinprod_jac21, &
       sll_f_sinprod_jac22, &
       (/0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/) )

  print *, 'initialization of the interpolators'
 ! Set up the interpolators for the field
  call interp_eta1%initialize( nc_eta1+1, 0.0_f64, 1.0_f64, sll_p_periodic )
  call interp_eta2%initialize( nc_eta2+1, 0.0_f64, 1.0_f64, sll_p_periodic )
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2


  print*, 'initialization of sll_m_distribution_function'

  call init_landau%initialize(m,sll_p_cell_centered_field,0.001_f64)
  p_init_f => init_landau

  print*, 'landau initialized'

  call sll_s_initialize_distribution_function_2d( &
       df, &
       1.0_f64, &
       1.0_f64, &
       name, &
       m, &
       sll_p_cell_centered_field, &
       interp_eta1_ptr, &
       interp_eta2_ptr, &
       p_init_f )

  print*, 'write mesh and distribution function'

  istep = 0
  call sll_s_int2string(istep,cstep)
  df%name = trim(name)//cstep
  
  call sll_s_write_scalar_field_2d(df,multiply_by_jacobian=.true.) 

!!$  x1_min =  0.0_f64; x1_max =  2.0_f64 * sll_p_pi
!!$  x2_min =  0.0_f64; x2_max =  2.0_f64 * sll_p_pi
!!$
!!$  v1_min = -6.0_f64; v1_max =  6.0_f64 
!!$  v2_min = -6.0_f64; v2_max =  6.0_f64 
!!$
!!$  nc_v1 = 31; nc_v2 = 31
!!$  nc_x1 = 31; nc_x2 = 31
!!$  
!!$  geom_x => new_geometry_2D('cartesian')
!!$  geom_v => new_geometry_2D('cartesian')
!!$
!!$  mesh_x => new_mesh_descriptor_2D(x1_min, x1_max, nc_x1, &
!!$          PERIODIC, x2_min, x2_max, nc_x2, PERIODIC, geom_x)
!!$
!!$  mesh_v => new_mesh_descriptor_2D(v1_min, v1_max, nc_v1, &
!!$          PERIODIC, v2_min, v2_max, nc_v2, PERIODIC, geom_v)
!!$  
!!$  call write_mesh_2D(mesh_x,"mesh_x")
!!$  call write_mesh_2D(mesh_v,"mesh_v")
!!$
!!$  dist_func_4D => sll_new_distribution_function_4D(mesh_x,mesh_v,NODE_CENTERED_DF,"df_4D")
!!$  call sll_init_distribution_function_4D( dist_func_4D, LANDAU)
!!$
!!$  call write_distribution_function (dist_func_4D) 
!!$
!!$  nnode_x1 = mesh_x%nc_eta1+1
!!$  nnode_v1 = mesh_v%nc_eta1+1
!!$  SLL_ALLOCATE(val(nnode_x1,nnode_v1), error)
!!$  do ix = 1, nnode_x1
!!$     do iv = 1, nnode_v1
!!$        val(ix,iv) = sum(dist_func_4D%field%data(ix,:,iv,:))
!!$     end do
!!$  end do
!!$
!!$  geom_xv => new_geometry_2D('cartesian')
!!$
!!$  mesh_xv => new_mesh_descriptor_2D(x1_min, x1_max, nc_x1, &
!!$             PERIODIC, v1_min, v1_max, nc_v1, PERIODIC, geom_xv)
!!$
!!$  call write_mesh_2D(mesh_xv,"mesh_xv")
!!$    
!!$  call write_vec1d(val,nc_x1+1,nc_v1+1,"df_on_xv_space","mesh_xv",NODE_CENTERED_DF)
    
  print *, 'Successful, exiting program.'
  
end program unit_test
