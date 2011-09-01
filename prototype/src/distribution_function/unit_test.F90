program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_mesh_types.h"
  use numeric_constants
  use distribution_function
  use sll_diagnostics
  implicit none
  
  sll_int32 :: nc_eta1, nc_eta2
  sll_real64 :: eta1_min, eta1_max, delta_eta1, delta_eta2, eta2_min, eta2_max
  type(geometry_2D), pointer :: geom
  type(mesh_descriptor_2D), pointer :: m2D_descriptor
  type(sll_distribution_function_2D_t), pointer :: dist_func
  character(32)  :: name = 'dist_func'

  nc_eta1 = 100
  nc_eta2 = 100

  print*, 'checking initialization of distribution_function'
  
  eta1_min = -6.0_f64
  eta1_max = 6.0_f64
  eta2_min = -6.0_f64
  eta2_max = 6.0_f64 

  geom => new_geometry_2D ('sinprod')
  m2D_descriptor => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, COMPACT, geom)
  dist_func => sll_new_distribution_function_2D(m2D_descriptor,name)
  delta_eta1 = get_df_delta_eta1 ( dist_func )
  delta_eta2 = get_df_delta_eta2 ( dist_func )
  call sll_init_distribution_function_2D( dist_func, GAUSSIAN, 'cell' )
  call write_mesh_2D(m2D_descriptor)
  call write_distribution_function (dist_func) 

  print *, 'Successful, exiting program.'
  
end program unit_test
