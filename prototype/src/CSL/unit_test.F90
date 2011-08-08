program unit_test
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
  use numeric_constants
  use distribution_function
  use advection_field
  use sll_diagnostics
  use sll_csl
  implicit none
  
  sll_int32 :: nc_eta1, nc_eta2, i1, i2, it
  sll_real64 :: eta1_min, eta1_max, delta_eta1, delta_eta2, eta2_min, eta2_max, x, v, deltat, val 
  type(geometry_2D), pointer :: geom
  type(mesh_descriptor_2D), pointer :: m2D_descriptor
  type(sll_distribution_function_2D_t), pointer :: dist_func
  type(field_2D_vec2), pointer :: rotating_field
  type(csl_workspace), pointer :: csl_work
  character(32),parameter  :: name = 'ions'

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

  dist_func => sll_new_distribution_function_2D(m2D_descriptor,name)
  
  delta_eta1 = get_df_delta_eta1 ( dist_func )
  delta_eta2 = get_df_delta_eta2 ( dist_func )
  call sll_init_distribution_function_2D( dist_func, GAUSSIAN )
  call write_mesh_2D(m2D_descriptor)
  call write_distribution_function ( dist_func )

  ! define rotating field
  rotating_field => new_field_2D_vec2(m2D_descriptor)
  x = eta1_min 
  do i1 = 1, nc_eta1+1
     v = eta2_min 
     do i2 = 1, nc_eta2+1
        FIELD_2D_AT_I_V1( rotating_field, i1, i2 ) = 1.0_f64
        FIELD_2D_AT_I_V2( rotating_field, i1, i2 ) = 0.0_f64
        v = v + delta_eta2
     end do
     x = x + delta_eta1
  end do
  ! initialize CSL  
  csl_work => new_csl_workspace( dist_func )
  ! run CSL method
  deltat = 0.1_f64
  do it = 1, 100
     call csl_first_order(csl_work, dist_func, rotating_field, deltat)
     call write_distribution_function ( dist_func )
  end do   
!!$  do i1 =1, nc_eta1 + 1
!!$     do i2 = 1, nc_eta2 +1
!!$        val = sll_get_df_val(dist_func, i1, i2)
!!$ !       print*, i1,i2, val
!!$     end do
!!$  end do
  print *, 'Successful, exiting program.'
  
end program unit_test
