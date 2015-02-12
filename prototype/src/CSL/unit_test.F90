program unit_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_field_2d.h"
#include "sll_memory.h"
#include "sll_constants.h"
#include "sll_file_io.h"

  use distribution_function
  use sll_csl
  !use sll_module_mapped_meshes_2d_cartesian
  use sll_common_coordinate_transformations
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_gaussian_2d_initializer
  use sll_cubic_spline_interpolator_1d
  implicit none
  
  sll_int32 :: nc_eta1_coarse, nc_eta2_coarse
  sll_int32 :: nc_eta1_fine!, nc_eta2_fine
  !sll_real64 :: delta_eta1_coarse, delta_eta2_coarse
  !sll_real64 :: delta_eta1_fine, delta_eta2_fine
  sll_real64 :: eta1_min, eta1_max,  eta2_min, eta2_max
  procedure(scalar_function_2D), pointer :: x1_coarse, x1_fine, x2_coarse, x2_fine, jac_coarse, jac_fine
  !class(sll_mapped_mesh_2d_base), pointer   :: m
  type(sll_logical_mesh_2d), pointer :: mesh_c, mesh_f
  class(sll_coordinate_transformation_2d_base), pointer   :: m
  class(scalar_field_2d_initializer_base), pointer    :: p_init_f
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
  !type(sll_mapped_mesh_2d_cartesian),target :: mesh_c,mesh_f
  !type(geometry_2D), pointer :: geomc, geomf
  !type(mesh_descriptor_2D), pointer :: coarse_mesh
  !type(mesh_descriptor_2D), pointer :: fine_mesh
  type(sll_distribution_function_2D) :: dist_func_fine
  type(sll_distribution_function_2D) :: dist_func_coarse
  !type(scalar_field_2D), pointer :: rotating_field
  !type(scalar_field_2D), pointer :: uniform_field
  !type(csl_workspace), pointer :: csl_work
  character(32),parameter  :: name = 'distribution_function'
  type(init_gaussian_2d),target :: pgaussian
  type(cubic_spline_1d_interpolator), target  :: interp_eta1
  type(cubic_spline_1d_interpolator), target  :: interp_eta2

  
  eta1_min =  -8.0_f64
  eta1_max =  8.0_f64
  eta2_min =  -8.0_f64
  eta2_max =  8.0_f64 !2*sll_pi ! 8.0_f64
  nc_eta1_coarse = 100
  nc_eta2_coarse = 100
  
  
  mesh_c => new_logical_mesh_2d( &
       nc_eta1_coarse, &
       nc_eta2_coarse,  &
       eta1_min,       &
       eta1_max,       &
       eta2_min,       &
       eta2_max       &
   )
  m => new_coordinate_transformation_2d_analytic( &
       "mesh2d_cart",      &
       mesh_c,             &
       identity_x1,    &
       identity_x2,    &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       (/0.0_f64/) ) 

  call pgaussian%initialize( m, CELL_CENTERED_FIELD)

  p_init_f => pgaussian
  
  ! Set up the interpolators for the distribution function
  call interp_eta1%initialize( nc_eta1_coarse+1, 0.0_f64, 1.0_f64, SLL_PERIODIC )
  call interp_eta2%initialize( nc_eta2_coarse+1, 0.0_f64, 1.0_f64, SLL_PERIODIC )
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2

  !geomc => new_geometry_2D ('cartesian', nc_eta1_coarse,  nc_eta2_coarse)
  !nc_eta1_fine = 200
  !nc_eta2_fine = 200
  !geomf => new_geometry_2D ('cartesian', nc_eta1_fine,  nc_eta2_fine)
  
  !coarse_mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1_coarse, &
  !     COMPACT, eta2_min, eta2_max, nc_eta2_coarse, COMPACT, geomc)
  !fine_mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1_fine, &
  !     PERIODIC, eta2_min, eta2_max, nc_eta2_fine, PERIODIC, geomf)
  !dist_func_coarse => sll_new_distribution_function_2D(coarse_mesh,CELL_CENTERED_DF, name)
  !dist_func_fine => sll_new_distribution_function_2D(fine_mesh,CELL_CENTERED_DF,name)

  !delta_eta1_coarse = get_df_delta_eta1 ( dist_func_coarse )
  !delta_eta2_coarse = get_df_delta_eta2 ( dist_func_coarse )
  !delta_eta1_fine   = get_df_delta_eta1 ( dist_func_fine )
  !delta_eta2_fine   = get_df_delta_eta2 ( dist_func_fine )

  !x1_coarse => get_df_x1 ( dist_func_coarse )
  !x2_coarse => get_df_x2 ( dist_func_coarse )
  !jac_coarse => get_df_jac ( dist_func_coarse )
  !x1_fine => get_df_x1 ( dist_func_fine )
  !x2_fine => get_df_x2 ( dist_func_fine )
  !jac_fine => get_df_jac ( dist_func_fine )


  !call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )

  call initialize_distribution_function_2d( &
    dist_func_coarse, &
    1._f64, &
    1._f64, &
    "f_coarse", &
    m, &
    CELL_CENTERED_FIELD, &
    interp_eta1_ptr, &
    interp_eta2_ptr, &
    p_init_f )



  !call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
  !call write_mesh_2D(fine_mesh)

  !call write_distribution_function ( dist_func_fine )

  !call sll_init_distribution_function_2D( dist_func_coarse, GAUSSIAN )

end program unit_test 

!  Print*, 'checking advection of a Gaussian in a uniform field'
!  !    no splitting error. First and second order splitting should be same.
!  !    only interpolation error
!  
!  ! define uniform field on coarse_mesh (using stream function)
!  uniform_field => new_field_2D_vec1(coarse_mesh)
!  ! components of field
!  alpha1 = -1.0_f64
!  alpha2 = -1.0_f64
!  eta1 = eta1_min 
!  do i1 = 1, nc_eta1_coarse+1
!     eta2 = eta2_min 
!     do i2 = 1, nc_eta2_coarse+1
!        !FIELD_2D_AT_I_V1( uniform_field, i1, i2 ) = -1.0_f64
!        !FIELD_2D_AT_I_V2( uniform_field, i1, i2 ) = -1.0_f64
!        FIELD_2D_AT_I( uniform_field, i1, i2 ) = alpha1 * x2_coarse(eta1,eta2) - alpha2 * x1_coarse(eta1,eta2)
!        eta2 = eta2 + delta_eta2_coarse
!     end do
!     eta1 = eta1 + delta_eta1_coarse
!  end do
!  ! initialize CSL  
!  csl_work => new_csl_workspace( dist_func_coarse )
!  ! run CSL method for 10 time steps
!  n_steps = 10
!  deltat = 1.0_f64/n_steps
!    do it = 1, n_steps
!     call csl_first_order(csl_work, dist_func_coarse, uniform_field, deltat)
!     !call csl_second_order(csl_work, dist_func, rotating_field, rotating_field, deltat)
!     !call write_distribution_function ( dist_func_coarse )
!  end do   
!  ! compute error when Gaussian arrives at center (t=1)
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_coarse  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_coarse
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_coarse  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_coarse 
!        val = sll_get_df_val(dist_func_coarse, i1, i2) / jac_coarse(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*(x1_coarse(eta1,eta2)**2+x2_coarse(eta1,eta2)**2))))
!        eta2 = eta2 + delta_eta2_coarse
!     end do
!     eta1 = eta1 + delta_eta1_coarse
!  end do
!  print*, '    coarse mesh, 1st order splitting, 100 cells, 10 time steps. Error= ', error
!!!$
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_coarse, GAUSSIAN)
!  ! run CSL method using 20 time steps
!  n_steps = 20
!  deltat = 1.0_f64/n_steps
!    do it = 1, n_steps
!     call csl_first_order(csl_work, dist_func_coarse, uniform_field, deltat)
!     !call csl_second_order(csl_work, dist_func_coarse, uniform_field, uniform_field, deltat)
!     !call write_distribution_function ( dist_func )
!  end do   
!  ! compute error when Gaussian arrives at center (t=1)
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_coarse  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_coarse
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_coarse  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_coarse 
!        val = sll_get_df_val(dist_func_coarse, i1, i2) / jac_coarse(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*(x1_coarse(eta1,eta2)**2+x2_coarse(eta1,eta2)**2))))
!        eta2 = eta2 + delta_eta2_coarse
!     end do
!     eta1 = eta1 + delta_eta1_coarse
!  end do
!  print*, '    coarse mesh, 1st order splitting, 100 cells, 20 time steps. Error= ', error
!  
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_coarse, GAUSSIAN)
!  ! run CSL method using 10 time steps and second order splitting
!  n_steps = 10
!  deltat = 1.0_f64/n_steps
!    do it = 1, n_steps
!       call csl_second_order(csl_work, dist_func_coarse, uniform_field, uniform_field, deltat)
!    end do
!  ! compute error when Gaussian arrives at center (t=1)
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_coarse  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_coarse
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_coarse  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_coarse 
!        val = sll_get_df_val(dist_func_coarse, i1, i2) / jac_coarse(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*(x1_coarse(eta1,eta2)**2+x2_coarse(eta1,eta2)**2))))
!        eta2 = eta2 + delta_eta2_coarse
!     end do
!     eta1 = eta1 + delta_eta1_coarse
!  end do
!  print*, '    coarse mesh, 2nd order splitting, 100 cells, 10 time steps. Error= ', error
!
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_coarse, GAUSSIAN )
!  ! run CSL method using 20 time steps and second order splitting
!  n_steps = 20
!  deltat = 1.0_f64/n_steps
!    do it = 1, n_steps
!       call csl_second_order(csl_work, dist_func_coarse, uniform_field, uniform_field, deltat)
!    end do
!  ! compute error when Gaussian arrives at center (t=1)
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_coarse  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_coarse
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_coarse  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_coarse 
!        val = sll_get_df_val(dist_func_coarse, i1, i2) / jac_coarse(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*(x1_coarse(eta1,eta2)**2+x2_coarse(eta1,eta2)**2))))
!        eta2 = eta2 + delta_eta2_coarse
!     end do
!     eta1 = eta1 + delta_eta1_coarse
!  end do
!  print*, '    coarse mesh, 2nd order splitting, 100 cells, 20 time steps. Error= ', error
!
!  ! define uniform field on fine_mesh
!  call delete_field_2D_vec1(uniform_field)
!  uniform_field => new_field_2D_vec1(fine_mesh)
!  alpha1 = -1.0_f64
!  alpha2 = -1.0_f64
!  eta1 = eta1_min 
!  do i1 = 1, nc_eta1_fine+1
!     eta2 = eta2_min 
!     do i2 = 1, nc_eta2_fine+1
!        !FIELD_2D_AT_I_V1( uniform_field, i1, i2 ) = -1.0_f64
!        !FIELD_2D_AT_I_V2( uniform_field, i1, i2 ) = -0.0_f64
!        FIELD_2D_AT_I( uniform_field, i1, i2 ) = alpha1 * x2_fine(eta1,eta2) - alpha2 * x1_fine(eta1,eta2)
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  ! reinitialize CSL  
!  !call delete_csl_workspace( csl_work )
!  csl_work => new_csl_workspace( dist_func_fine )
!  print*, 'working now on fine mesh'
!
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
!  ! run CSL method for 10 time steps
!  n_steps = 10
!  deltat = 1.0_f64/n_steps
!  !deltat = 0.1_f64 * delta_eta1_fine
!  !print*, 'deltat=', deltat, 'delta_eta1=',delta_eta1_fine, 'delta_eta2=',delta_eta2_fine
!    do it = 1, n_steps
!       call csl_first_order(csl_work, dist_func_fine, uniform_field, deltat)
!       !call csl_second_order(csl_work, dist_func, rotating_field, rotating_field, deltat)
!       !call write_distribution_function ( dist_func_fine )
!  end do   
!  ! compute error when Gaussian arrives at center (t=1)
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_fine  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_fine
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_fine  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_fine 
!        val = sll_get_df_val(dist_func_fine, i1, i2) / jac_fine(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*(x1_fine(eta1,eta2)**2+x2_fine(eta1,eta2)**2))))
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  print*, '    fine mesh,   1st order splitting, 200 cells, 10 time steps. Error= ', error
!
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
!  ! run CSL method using 20 time steps
!  n_steps = 20
!  deltat = 1.0_f64/n_steps
!    do it = 1, n_steps
!     call csl_first_order(csl_work, dist_func_fine, uniform_field, deltat)
!     !call csl_second_order(csl_work, dist_func, rotating_field, rotating_field, deltat)
!     !call write_distribution_function ( dist_func )
!  end do   
!  ! compute error when Gaussian arrives at center (t=1)
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_fine  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_fine
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_fine  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_fine 
!        val = sll_get_df_val(dist_func_fine, i1, i2) / jac_fine(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*(x1_fine(eta1,eta2)**2+x2_fine(eta1,eta2)**2))))
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  print*, '    fine mesh,   1st order splitting, 200 cells, 20 time steps. Error= ', error
!  
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
!  ! run CSL method using 10 time steps and second order splitting
!  n_steps = 10
!  deltat = 1.0_f64/n_steps
!    do it = 1, n_steps
!       call csl_second_order(csl_work, dist_func_fine, uniform_field, uniform_field, deltat)
!    end do
!  ! compute error when Gaussian arrives at center (t=1)
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_fine  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_fine
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_fine  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_fine 
!        val = sll_get_df_val(dist_func_fine, i1, i2) / jac_fine(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*(x1_fine(eta1,eta2)**2+x2_fine(eta1,eta2)**2))))
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  print*, '    fine mesh,   2nd order splitting, 200 cells, 10 time steps. Error= ', error
!
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
!  ! run CSL method using 20 time steps and second order splitting
!  n_steps = 20
!  deltat = 1.0_f64/n_steps
!    do it = 1, n_steps
!       call csl_second_order(csl_work, dist_func_fine, uniform_field, uniform_field, deltat)
!    end do
!  ! compute error when Gaussian arrives at center (t=1)
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_fine  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_fine
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_fine  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_fine 
!        val = sll_get_df_val(dist_func_fine, i1, i2) / jac_fine(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*(x1_fine(eta1,eta2)**2+x2_fine(eta1,eta2)**2))))
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  print*, '    fine mesh,   2nd order splitting, 200 cells, 20 time steps. Error= ', error
!  print*, '    Conclusions: no splitting error, no time integration error, only interpolation error (4th order)'
!
!  print*, 'checking advection in rotating field' 
!  ! define rotating field
!  rotating_field => new_field_2D_vec1(fine_mesh)
!  eta1 = eta1_min 
!  do i1 = 1, nc_eta1_fine+1
!     eta2 = eta2_min 
!     do i2 = 1, nc_eta2_fine+1
!        !FIELD_2D_AT_I_V1( rotating_field, i1, i2 ) = x2_fine(eta1,eta2) 
!        !FIELD_2D_AT_I_V2( rotating_field, i1, i2 ) = -x1_fine(eta1,eta2) 
!        FIELD_2D_AT_I( rotating_field, i1, i2 ) = 0.5_f64*(x1_coarse(eta1,eta2)**2 + x2_coarse(eta1,eta2)**2)
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
!  ! run CSL method
!  n_steps = 10
!  deltat = 0.5_f64*sll_pi/n_steps  ! do one quarter turn
!  do it = 1, n_steps
!     call csl_first_order(csl_work, dist_func_fine, rotating_field, deltat)
!  end do
!  ! compute error after one quarter turn
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_fine  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_fine
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_fine  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_fine 
!        val = sll_get_df_val(dist_func_fine, i1, i2) / jac_fine(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*((x1_fine(eta1,eta2)-1.0_f64)**2 &
!             + (x2_fine(eta1,eta2)+1.0_f64)**2))))
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  print*, '    fine mesh, 1st order splitting, 200 cells, 10 time steps,  Error=', error
!  error1 = error
!
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
!  ! run CSL method
!  n_steps = 20
!  deltat = 0.5_f64*sll_pi/n_steps  ! do one quarter turn
!  do it = 1, n_steps
!     call csl_first_order(csl_work, dist_func_fine, rotating_field, deltat)
!  end do
!  ! compute error after one quarter turn
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_fine  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_fine
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_fine  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_fine 
!        val = sll_get_df_val(dist_func_fine, i1, i2) / jac_fine(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*((x1_fine(eta1,eta2)-1.0_f64)**2 &
!             +(x2_fine(eta1,eta2)+1.0_f64)**2))))
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  print*, '    fine mesh, 1st order splitting, 200 cells, 20 time steps,  Error=', error
!  print*, '                   order=', error1/error , ' should be close to 2'
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
!  ! run CSL method
!  n_steps = 10
!  deltat = 0.5_f64*sll_pi/n_steps  ! do one quarter turn
!  do it = 1, n_steps
!     call csl_second_order(csl_work, dist_func_fine, rotating_field, rotating_field, deltat)
!     call write_distribution_function ( dist_func_fine )
!  end do
!  ! compute error after one turn
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_fine  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_fine
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_fine  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_fine 
!        val = sll_get_df_val(dist_func_fine, i1, i2) / jac_fine(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*((x1_fine(eta1,eta2)-1.0_f64)**2&
!             +(x2_fine(eta1,eta2)+1.0_f64)**2))))
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  print*, '    fine mesh, 2nd order splitting, 200 cells, 10 time steps,  Error=', error
!  error1 = error
!  ! reinitialize distribution function
!  call sll_init_distribution_function_2D( dist_func_fine, GAUSSIAN )
!  ! run CSL method
!  n_steps = 20
!  deltat = 0.5_f64*sll_pi/n_steps  ! do one quarter turn
!  do it = 1, n_steps
!     call csl_second_order(csl_work, dist_func_fine, rotating_field, rotating_field, deltat)
!  end do
!  ! compute error after one turn
!  error = 0.0_f64
!  eta1 = eta1_min + 0.5_f64 * delta_eta1_fine  ! eta1 at midpoint of cell
!  do i1 = 1, nc_eta1_fine
!     eta2 = eta2_min + 0.5_f64 * delta_eta2_fine  ! eta2 at midpoint of cell
!     do i2 = 1, nc_eta2_fine 
!        val = sll_get_df_val(dist_func_fine, i1, i2) / jac_fine(eta1,eta2)
!        error = max (error, abs(val - exp(-0.5_f64*((x1_fine(eta1,eta2)-1.0_f64)**2 &
!             +(x2_fine(eta1,eta2)+1.0_f64)**2))))
!        eta2 = eta2 + delta_eta2_fine
!     end do
!     eta1 = eta1 + delta_eta1_fine
!  end do
!  print*, '    fine mesh, 2nd order splitting, 200 cells, 20 time steps,  Error=', error
!  print*, '                   order=', error1/error, ' should be close to 4'
!  print*, '    Conclusion. Time splitting error of right order dominates. '
!  print*, '                No time integration error on 1D ode solvers'
!
!
!  ! might be good to implement another test case where each split step is not a constant coefficient advection
!  ! to also check the ode solver
!  print *, 'Successful, exiting program.'
!  
!end program unit_test
