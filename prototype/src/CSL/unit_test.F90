program unit_test
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
  use numeric_constants
  use distribution_function
  use advection_field
  use sll_diagnostics
  implicit none
  
  sll_int32 :: nc_eta1, nc_eta2, i1, i2
  sll_real64 :: eta1_min, eta1_max, delta_eta1, delta_eta2, eta2_min, eta2_max, x, v
  type(geometry_2D), pointer :: geom
  type(mesh_descriptor_2D), pointer :: m2D_descriptor
  type(sll_distribution_function_2D_t), pointer :: dist_func
  type(field_2D_vec2), pointer :: rotating_field

  print*, 'checking advection in rotating field'
  nc_eta1 = 64
  nc_eta2 = 64
  eta1_min = -6.0_f64
  eta1_max =  6.0_f64
  eta2_min = -6.0_f64
  eta2_max =  6.0_f64

  geom => new_geometry_2D ()
  m2D_descriptor => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, COMPACT, geom)
  dist_func => sll_new_distribution_function_2D(m2D_descriptor)
  delta_eta1 = get_df_delta_eta1 ( dist_func )
  delta_eta2 = get_df_delta_eta2 ( dist_func )
  call sll_init_distribution_function_2D( dist_func, GAUSSIAN )
  call write_distribution_function ( dist_func )
  
  ! define rotating field
  rotating_field => new_field_2D_vec2(m2D_descriptor)
  x = eta1_min 
  do i1 = 1, nc_eta1+1
     v = eta2_min 
     do i2 = 1, nc_eta2+1
        FIELD_2D_AT_I_V1( rotating_field, i1, i2 ) = v
        FIELD_2D_AT_I_V2( rotating_field, i1, i2 ) = -x
        v = v + delta_eta2
     end do
     x = x + delta_eta1
  end do
  print *, 'Successful, exiting program.'
  
end program unit_test
