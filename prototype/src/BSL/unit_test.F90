program unit_test
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"

  use numeric_constants
  use distribution_function
  use sll_diagnostics
  use sll_bsl

  implicit none


  call test_1d()
  call test_2d()
  
  contains

  subroutine test_1d()

  sll_int32  :: nc_eta1, nc_eta2
  sll_real64 :: delta_eta1, delta_eta2
  sll_int32  :: i1, i2, it, n_steps
  sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
  sll_real64 :: eta1, eta2, delta_t, val, error

  type(mesh_descriptor_1D),             pointer :: mesh_x
  type(mesh_descriptor_1D),             pointer :: mesh_v
  type(mesh_descriptor_2D),             pointer :: mesh_xv
  type(sll_distribution_function_2D_t), pointer :: df_2d
  type(field_1D_vec1),                  pointer :: uniform_field_x
  type(field_1D_vec1),                  pointer :: uniform_field_v
  type(bsl_workspace_1d),               pointer :: bsl_work_x
  type(bsl_workspace_1d),               pointer :: bsl_work_v

  print*,'*********************'
  print*,' 1D case             '
  print*,' 1D in x and 1D in v '
  print*,'*********************'

  print*, 'set domain size'
  eta1_min =  -8.0_f64; eta1_max =  8.0_f64
  eta2_min =  -8.0_f64; eta2_max =  8.0_f64 

  nc_eta1 = 100; nc_eta2 = 100

  print*, 'create 1d meshes in x and v'
  mesh_x => new_mesh_descriptor_1D(eta1_min, eta1_max, nc_eta1, PERIODIC)
  mesh_v => new_mesh_descriptor_1D(eta2_min, eta2_max, nc_eta2, PERIODIC)

  print*, 'create 2d mesh from mesh_x and mesh_v'
  mesh_xv => mesh_x * mesh_v

  print*, 'create 2d distribution function f(x,v)'
  df_2d => sll_new_distribution_function_2D(mesh_xv, NODE_CENTERED_DF, 'one_d')

  delta_eta1 = mesh_x%delta_eta1
  delta_eta2 = mesh_v%delta_eta1

  call write_mesh_2D(mesh_xv)

  print*, 'initialize 2d distribution function f(x,v) gaussian'
  call sll_init_distribution_function_2D( df_2d, GAUSSIAN)
  call write_distribution_function ( df_2d )

  Print*, 'checking advection of a Gaussian in a uniform field'
  
  uniform_field_x => new_field_1D_vec1(mesh_x)
  uniform_field_v => new_field_1D_vec1(mesh_v)

  eta1 = eta1_min 
  do i1 = 1, nc_eta1+1
     uniform_field_x%data(i1) = 1_f64 
     eta1 = eta1 + delta_eta1
  end do

  eta2 = eta2_min 
  do i2 = 1, nc_eta2+1
     uniform_field_v%data(i2) = 1_f64 
     eta2 = eta2 + delta_eta2
  end do
  
  bsl_work_x => new_bsl_workspace(uniform_field_x)
  bsl_work_v => new_bsl_workspace(uniform_field_v)

  ! run BSL method using 10 time steps and second order splitting
  n_steps = 10
  delta_t = 1.0_f64/n_steps
  do it = 1, n_steps
     do i2 = 1, nc_eta2
        call bsl_step_1d( bsl_work_x, df_2d%field%data(:,i2), &
                          uniform_field_x, delta_t )
     end do
     do i1 = 1, nc_eta1
        call bsl_step_1d( bsl_work_v, df_2d%field%data(i1,:), &
                          uniform_field_v, delta_t )
     end do

     call write_distribution_function ( df_2d )

  end do
    
  ! compute error when Gaussian arrives at center (t=1)

  error = 0.0_f64
  eta1 = eta1_min
  do i1 = 1, nc_eta1
     eta2 = eta2_min
     do i2 = 1, nc_eta2
        val = sll_get_df_val(df_2d, i1, i2) 
        error = max(error,abs(val-exp(-0.5_f64*(eta1*eta1+eta2*eta2))))
        eta2 = eta2 + delta_eta2
     end do
     eta1 = eta1 + delta_eta1
  end do

  print*, ' 100 nodes, 10 time steps. Error= ', error

  print *, 'Successful, exiting program.'

  end subroutine test_1d
  
  subroutine test_2d()

  sll_int32  :: nc_eta1, nc_eta2
  sll_real64 :: delta_eta1, delta_eta2
  sll_int32  :: i1, i2, it, n_steps
  sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
  sll_real64 :: eta1, eta2, delta_t, val, error

  type(geometry_2D),                    pointer :: geom
  type(mesh_descriptor_2D),             pointer :: mesh
  type(sll_distribution_function_2D_t), pointer :: dist_func
  type(field_2D_vec1),                  pointer :: uniform_field
  type(bsl_workspace_2d),               pointer :: bsl_work

  print*,'*********'
  print*,' 2D case '
  print*,'*********'
  eta1_min =  -8.0_f64; eta1_max =  8.0_f64
  eta2_min =  -8.0_f64; eta2_max =  8.0_f64 

  geom => new_geometry_2D ('cartesian')

  nc_eta1 = 100; nc_eta2 = 100

  mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
          PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

  dist_func => sll_new_distribution_function_2D(mesh, NODE_CENTERED_DF, 'two_d')

  delta_eta1 = get_df_delta_eta1 (dist_func)
  delta_eta2 = get_df_delta_eta2 (dist_func)

  call write_mesh_2D(mesh)
  call sll_init_distribution_function_2D( dist_func, GAUSSIAN)
  call write_distribution_function ( dist_func )

  Print*, 'checking advection of a Gaussian in a uniform field'
  !    no splitting error. First and second order splitting should be same.
  !    only interpolation error
  
  ! define uniform field on mesh (using stream function)
  uniform_field => new_field_2D_vec1(mesh)
  ! components of field
  eta1 = eta1_min 
  do i1 = 1, nc_eta1+1
     eta2 = eta2_min 
     do i2 = 1, nc_eta2+1
        FIELD_2D_AT_I( uniform_field, i1, i2 ) = 1_f64 
        eta2 = eta2 + delta_eta2
     end do
     eta1 = eta1 + delta_eta1
  end do
  ! initialize BSL  
  bsl_work => new_bsl_workspace( dist_func)

  ! run BSL method using 10 time steps and second order splitting
  n_steps = 10
  delta_t = 1.0_f64/n_steps
  do it = 1, n_steps
     call bsl_step_2d( bsl_work, dist_func, uniform_field, delta_t )
     call write_distribution_function ( dist_func )
  end do
    
  ! compute error when Gaussian arrives at center (t=1)

  error = 0.0_f64
  eta1 = eta1_min
  do i1 = 1, nc_eta1
     eta2 = eta2_min
     do i2 = 1, nc_eta2
        val = sll_get_df_val(dist_func, i1, i2) 
        error = max(error,abs(val-exp(-0.5_f64*(eta1*eta1+eta2*eta2))))
        eta2 = eta2 + delta_eta2
     end do
     eta1 = eta1 + delta_eta1
  end do

  print*, ' 2nd order splitting, 100 nodes, 10 time steps. Error= ', error

  ! might be good to implement another test case where each 
  ! split step is not a constant coefficient advection
  ! to also check the ode solver
  print *, 'Successful, exiting program.'

  end subroutine test_2d
  
end program unit_test
