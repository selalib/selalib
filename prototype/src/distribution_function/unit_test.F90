program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_mesh_2d.h"
  use numeric_constants
  use distribution_function
  use sll_diagnostics
  use geometry_functions
  use initial_distribution_functions
  implicit none
 
  sll_int32 :: nc_eta1, nc_eta2
  type(mesh_2d_analytic) :: mesh2d
  type(sll_distribution_function_2D_t) :: dist_func
  character(32)  :: name = 'dist_func'
  procedure(scalar_function_2D), pointer :: p_init_f, px1, px2, pjac
  sll_int32  :: ierr
  sll_int32 :: ix, iv, nnode_x1, nnode_v1


  nc_eta1 = 100
  nc_eta2 = 100

  print*, 'initialization of mesh'
  
  px1 => sinprod_x1
  px2 => sinprod_x2
  pjac => sinprod_jac
  call new_mesh_2d_analytic ( mesh2d, nc_eta1, nc_eta2, px1, px2, pjac)

  print*, 'initialization of distribution_function'

  p_init_f => gaussian
  call sll_new_distribution_function_2D(dist_func,mesh2d,CELL_CENTER_FIELD, &
       name, p_init_f)
 
  print*, 'write mesh and distribution function'
  call write_mesh_2d(mesh2d)
  call write_scalar_field_2d(dist_func,name,jacobian=.true.) 

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
