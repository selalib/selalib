program unit_test
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"

  use numeric_constants
  use distribution_function
  use sll_diagnostics
  use sll_bsl

  implicit none
  
  sll_int32  :: nc_eta1, nc_eta2
  sll_real64 :: delta_eta1, delta_eta2
  sll_int32  :: i1, i2, it, n_steps
  sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
  sll_real64 :: eta1, eta2, deltat, val, error
  sll_real64 :: alpha1, alpha2

  procedure(scalar_function_2D), pointer :: x1, x2

  type(geometry_2D), pointer :: geom
  type(mesh_descriptor_2D), pointer :: mesh

  type(sll_distribution_function_2D_t), pointer :: dist_func

  type(field_2D_vec1), pointer :: uniform_field
  type(bsl_workspace), pointer :: bsl_work

  eta1_min =  -8.0_f64; eta1_max =  8.0_f64
  eta2_min =  -8.0_f64; eta2_max =  8.0_f64 

  geom => new_geometry_2D ('cartesian')

  nc_eta1 = 100; nc_eta2 = 100

  mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
          PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

  dist_func => sll_new_distribution_function_2D(mesh, NODE_CENTERED_DF, 'f')

  delta_eta1 = get_df_delta_eta1 (dist_func)
  delta_eta2 = get_df_delta_eta2 (dist_func)

  x1 => get_df_x1 ( dist_func )
  x2 => get_df_x2 ( dist_func )
  
  call write_mesh_2D(mesh)
  call sll_init_distribution_function_2D( dist_func, GAUSSIAN)
  call write_distribution_function ( dist_func )

  Print*, 'checking advection of a Gaussian in a uniform field'
  !    no splitting error. First and second order splitting should be same.
  !    only interpolation error
  
  ! define uniform field on mesh (using stream function)
  uniform_field => new_field_2D_vec1(mesh)
  ! components of field
  alpha1 = -1.0_f64
  alpha2 = -1.0_f64
  eta1 = eta1_min 
  do i1 = 1, nc_eta1+1
     eta2 = eta2_min 
     do i2 = 1, nc_eta2+1
        !FIELD_2D_AT_I_V1( uniform_field, i1, i2 ) = -1.0_f64
        !FIELD_2D_AT_I_V2( uniform_field, i1, i2 ) = -1.0_f64
        FIELD_2D_AT_I( uniform_field, i1, i2 ) = alpha1 * x2(eta1,eta2) - alpha2 * x1(eta1,eta2)
        eta2 = eta2 + delta_eta2
     end do
     eta1 = eta1 + delta_eta1
  end do
  ! initialize BSL  
  bsl_work => new_bsl_workspace( dist_func)

  ! run BSL method using 10 time steps and second order splitting
  n_steps = 10
  deltat = 1.0_f64/n_steps
  do it = 1, n_steps
     call bsl_second_order(bsl_work, dist_func, uniform_field, uniform_field, deltat)
     call write_distribution_function ( dist_func )
  end do
    
  ! compute error when Gaussian arrives at center (t=1)

  error = 0.0_f64
  eta1 = eta1_min
  do i1 = 1, nc_eta1
     eta2 = eta2_min
     do i2 = 1, nc_eta2
        val = sll_get_df_val(dist_func, i1, i2) 
        error = max (error, abs(val - exp(-0.5_f64*(x1(eta1,eta2)**2+x2(eta1,eta2)**2))))
        eta2 = eta2 + delta_eta2
     end do
     eta1 = eta1 + delta_eta1
  end do

  print*, ' 2nd order splitting, 100 nodes, 10 time steps. Error= ', error

  ! might be good to implement another test case where each split step is not a constant coefficient advection
  ! to also check the ode solver
  print *, 'Successful, exiting program.'
  
end program unit_test
