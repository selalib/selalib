program unit_test

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
  use numeric_constants
  use distribution_function
  use sll_module_mapped_meshes_2d
  use user_geometry_functions
  use geometry_functions
  use initial_distribution_functions
  implicit none

  sll_int32 :: nc_eta1, nc_eta2
  type(sll_mapped_mesh_2d_analytic) :: mesh2d
  type(sll_distribution_function_2d)   :: df 
  character(32)  :: name = 'dist_func'
  character(len=4) :: cstep
  procedure(scalar_function_2D), pointer :: p_init_f 
  sll_int32  :: ierr, istep
  sll_int32 :: ix, iv, nnode_x1, nnode_v1

  nc_eta1 = 100
  nc_eta2 = 100

  print*, 'initialization of mesh'
  ! functions defining the mesh are defined in user_geometry_functions.F90
  call mesh2d%initialize( &
       "mesh2d",  &
       nc_eta1+1, &
       nc_eta2+1, &
       x1_cartesian, &
       x2_cartesian, &
       jac1_cartesian, &
       zero_function, &
       zero_function, &
       jac2_cartesian )

print*, 'initialization of distribution_function'

  p_init_f => gaussian
!  call sll_new_distribution_function_2D(df,mesh2d,NODE_CENTERED_FIELD, &
!       name, p_init_f)

  call initialize_distribution_function_2d( &
       df, &
       1.0_f64, &
       1.0_f64, &
       name, &
       mesh2d, &
       NODE_CENTERED_FIELD, &
       p_init_f )
  print*, 'write mesh and distribution function'

  istep = 0
  call int2string(istep,cstep)
  df%name = trim(name)//cstep
  call write_scalar_field_2d(df,multiply_by_jacobian=.true.) 

  print *, 'Successful, exiting program.'

end program unit_test
