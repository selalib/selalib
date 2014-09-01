program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
  use sll_constants
  use distribution_function
  use sll_common_coordinate_transformations
  !use sll_module_mapped_meshes_2d
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_landau_2d_initializer
  use sll_cubic_spline_interpolator_1d
  implicit none
 
  sll_int32 :: nc_eta1, nc_eta2
!  type(sll_mapped_mesh_2d_analytic), target :: mesh2d
  !class(sll_mapped_mesh_2d_base), pointer   :: m
  class(sll_coordinate_transformation_2d_base), pointer   :: m
  type(sll_logical_mesh_2d), pointer :: m_log
  class(scalar_field_2d_initializer_base), pointer    :: p_init_f
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
  type(sll_distribution_function_2d)   :: df 
  character(32)  :: name = 'dist_func'
  character(len=4) :: cstep
  type(init_landau_2d), target :: init_landau
  type(cubic_spline_1d_interpolator), target  :: interp_eta1
  type(cubic_spline_1d_interpolator), target  :: interp_eta2

  sll_int32 :: istep

  nc_eta1 = 100
  nc_eta2 = 100

  print*, 'initialization of mesh'
  
 m_log => new_logical_mesh_2d( &
       nc_eta1, &
       nc_eta2  &
   )
  m => new_coordinate_transformation_2d_analytic( &
       "mesh2d_coll",      &
       m_log,             &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/) )

  print *, 'initialization of the interpolators'
 ! Set up the interpolators for the field
  call interp_eta1%initialize( nc_eta1+1, 0.0_f64, 1.0_f64, SLL_PERIODIC )
  call interp_eta2%initialize( nc_eta2+1, 0.0_f64, 1.0_f64, SLL_PERIODIC )
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2


  print*, 'initialization of distribution_function'

  call init_landau%initialize(m,CELL_CENTERED_FIELD,0.001_f64)
  p_init_f => init_landau

  print*, 'landau initialized'

  call initialize_distribution_function_2d( &
       df, &
       1.0_f64, &
       1.0_f64, &
       name, &
       m, &
       CELL_CENTERED_FIELD, &
       interp_eta1_ptr, &
       interp_eta2_ptr, &
       p_init_f )

  print*, 'write mesh and distribution function'

  istep = 0
  call int2string(istep,cstep)
  df%name = trim(name)//cstep
  
  call write_scalar_field_2d(df,multiply_by_jacobian=.true.) 

!!$  x1_min =  0.0_f64; x1_max =  2.0_f64 * sll_pi
!!$  x2_min =  0.0_f64; x2_max =  2.0_f64 * sll_pi
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
