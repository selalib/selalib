program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_mesh_types.h"
  use numeric_constants
  use distribution_function
  use sll_diagnostics

  implicit none

  print*, 'checking initialization of distribution_function'

  call test_2d()

  call test_2d_split()

  call test_4d()

  print *, 'Successful, exiting program.'
  
  contains

  subroutine test_2d_split()

  type(mesh_descriptor_1d), pointer :: mesh_x
  type(mesh_descriptor_1d), pointer :: mesh_v
  type(mesh_descriptor_2d), pointer :: mesh_xv

  type(sll_distribution_function_2D_t), pointer :: df_xv

  sll_real64 :: x_min, x_max
  sll_int32  :: nc_x
  sll_real64 :: v_min, v_max
  sll_int32  :: nc_v

  x_min =  0.0_f64 ; x_max = 4.0_f64 * sll_pi
  v_min = -6.0_f64 ; v_max = 6.0_f64

  nc_x = 64
  nc_v = 64

  mesh_x => new_mesh_descriptor_1d(x_min,x_max,nc_x, PERIODIC)
  mesh_v => new_mesh_descriptor_1d(v_min,v_max,nc_v, PERIODIC)

  mesh_xv =>  mesh_x * mesh_v

  call mesh_x%dump()
  call mesh_v%dump()

  call write_mesh_2D(mesh_xv)

  df_xv => sll_new_distribution_function_2D(mesh_xv,CELL_CENTERED_DF,"df_xv")

  call sll_init_distribution_function_2D( df_xv, GAUSSIAN)

  call write_distribution_function (df_xv) 

  end subroutine test_2d_split
  
  subroutine test_2d()

  sll_int32 :: nc_eta1, nc_eta2
  sll_real64 :: eta1_min, eta1_max, delta_eta1, delta_eta2, eta2_min, eta2_max
  type(geometry_2D), pointer :: geom
  type(mesh_descriptor_2D), pointer :: m2D_descriptor
  type(sll_distribution_function_2D_t), pointer :: dist_func
  character(32)  :: name = 'dist_func'

  nc_eta1 = 100
  nc_eta2 = 100
  
  eta1_min = -6.0_f64
  eta1_max = 6.0_f64
  eta2_min = -6.0_f64
  eta2_max = 6.0_f64 

  geom => new_geometry_2D ('sinprod')
  m2D_descriptor => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, COMPACT, geom)
  dist_func => sll_new_distribution_function_2D(m2D_descriptor,CELL_CENTERED_DF,name)
  delta_eta1 = get_df_delta_eta1 ( dist_func )
  delta_eta2 = get_df_delta_eta2 ( dist_func )
  call sll_init_distribution_function_2D( dist_func, GAUSSIAN)
  call write_mesh_2D(m2D_descriptor)
  call write_distribution_function (dist_func) 

  end subroutine test_2d
  
  subroutine test_4d()

  sll_int32  :: error
  sll_int32  :: nc_x1, nc_x2, nc_v1, nc_v2
  sll_real64 :: x1_min, x1_max, x2_min, x2_max
  sll_real64 :: v1_min, v1_max, v2_min, v2_max
  type(geometry_2D), pointer :: geom_x
  type(geometry_2D), pointer :: geom_v
  type(geometry_2D), pointer :: geom_xv
  type(mesh_descriptor_2D), pointer :: mesh_x
  type(mesh_descriptor_2D), pointer :: mesh_v
  type(mesh_descriptor_2D), pointer :: mesh_xv
  sll_real64, dimension(:,:), allocatable :: val
  sll_int32 :: ix, iv, nnode_x1, nnode_v1

  type(sll_distribution_function_4D_t), pointer :: dist_func_4D

  x1_min =  0.0_f64; x1_max =  2.0_f64 * sll_pi
  x2_min =  0.0_f64; x2_max =  2.0_f64 * sll_pi

  v1_min = -6.0_f64; v1_max =  6.0_f64 
  v2_min = -6.0_f64; v2_max =  6.0_f64 

  nc_v1 = 31; nc_v2 = 31
  nc_x1 = 31; nc_x2 = 31
  
  geom_x => new_geometry_2D('cartesian')
  geom_v => new_geometry_2D('cartesian')

  mesh_x => new_mesh_descriptor_2D(x1_min, x1_max, nc_x1, &
          PERIODIC, x2_min, x2_max, nc_x2, PERIODIC, geom_x)

  mesh_v => new_mesh_descriptor_2D(v1_min, v1_max, nc_v1, &
          PERIODIC, v2_min, v2_max, nc_v2, PERIODIC, geom_v)
  
  call write_mesh_2D(mesh_x,"mesh_x")
  call write_mesh_2D(mesh_v,"mesh_v")

  dist_func_4D => sll_new_distribution_function_4D(mesh_x,mesh_v,NODE_CENTERED_DF,"df_4D")
  call sll_init_distribution_function_4D( dist_func_4D, LANDAU)

  call write_distribution_function (dist_func_4D) 

  nnode_x1 = mesh_x%nc_eta1+1
  nnode_v1 = mesh_v%nc_eta1+1
  SLL_ALLOCATE(val(nnode_x1,nnode_v1), error)
  do ix = 1, nnode_x1
     do iv = 1, nnode_v1
        val(ix,iv) = sum(dist_func_4D%field%data(ix,:,iv,:))
     end do
  end do

  geom_xv => new_geometry_2D('cartesian')

  mesh_xv => new_mesh_descriptor_2D(x1_min, x1_max, nc_x1, &
             PERIODIC, v1_min, v1_max, nc_v1, PERIODIC, geom_xv)

  call write_mesh_2D(mesh_xv,"mesh_xv")
    
  call write_vec1d(val,nc_x1+1,nc_v1+1,"df_on_xv_space","mesh_xv",NODE_CENTERED_DF)
    

  end subroutine test_4d
  
end program unit_test
