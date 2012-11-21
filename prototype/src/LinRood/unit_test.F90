program unit_test
#include "sll_working_precision.h"
#include "sll_field_2d.h"
#include "sll_memory.h"
  use numeric_constants
  use geometry_functions
  use distribution_function
  use sll_scalar_field_2d
  use sll_module_mapped_meshes_2d_base
  use sll_module_mapped_meshes_2d
  use sll_scalar_field_initializers_base
  use sll_gaussian_2d_initializer
  use sll_linrood
  use sll_cubic_spline_interpolator_1d
  implicit none

  sll_int32 :: nc_eta1, nc_eta2
  type(sll_mapped_mesh_2d_analytic), target :: mesh2d, mesh2d_cart
  class(sll_mapped_mesh_2d_base), pointer   :: m
  type(sll_distribution_function_2d)   :: df 
  character(32)  :: name = 'dist_func'
  character(len=4) :: cstep
  sll_real64 :: deltat
  sll_int32 :: iter, nbiter 
  type(init_gaussian_2d), target :: init_gaussian
  class(scalar_field_2d_initializer_base), pointer    :: p_init_f
  type(scalar_field_2d) :: uniform_field, rotating_field, incompressible_field 
  type(linrood_plan) :: linrood
  type(cubic_spline_1d_interpolator), target  :: interp_eta1
  type(cubic_spline_1d_interpolator), target  :: interp_eta2
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
  ! interpolators for scalar field
  type(cubic_spline_1d_interpolator), target  :: interp_eta1_sf 
  type(cubic_spline_1d_interpolator), target  :: interp_eta2_sf
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr_sf
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr_sf
  ! interpolators for rotating field
  type(cubic_spline_1d_interpolator), target  :: interp_eta1_rf
  type(cubic_spline_1d_interpolator), target  :: interp_eta2_rf
  class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr_rf
  class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr_rf

  sll_int32  :: ierr, istep
  sll_int32 :: i1, i2
  sll_real64 :: alpha1, alpha2

  ! Define mapped mesh
  nc_eta1 = 100
  nc_eta2 = 100

  call mesh2d%initialize( &
       "mesh2d",      &
       nc_eta1+1,     &
       nc_eta2+1,     &
       sinprod_x1,    &
       sinprod_x2,    &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22 )

 call mesh2d_cart%initialize( &
       "mesh2d_cart",      &
       nc_eta1+1,     &
       nc_eta2+1,     &
       identity_x1,    &
       identity_x2,    &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22 )

 ! m => mesh2d
  m => mesh2d_cart

  print*, 'initialization of distribution_function'

  call init_gaussian%initialize( m, CELL_CENTERED_FIELD, 0.4_f64, 0.4_f64, 0.1_f64, 0.1_f64 )
  p_init_f => init_gaussian

 ! Set up the interpolators for the distribution function
  call interp_eta1%initialize( nc_eta1+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  call interp_eta2%initialize( nc_eta2+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  interp_eta1_ptr => interp_eta1
  interp_eta2_ptr => interp_eta2

 ! Set up the interpolators for the scalar field
  call interp_eta1_sf%initialize( nc_eta1+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  call interp_eta2_sf%initialize( nc_eta2+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  interp_eta1_ptr_sf => interp_eta1_sf
  interp_eta2_ptr_sf => interp_eta2_sf

 ! Set up the interpolators for the rotating field
  call interp_eta1_rf%initialize( nc_eta1+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  call interp_eta2_rf%initialize( nc_eta2+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )
  interp_eta1_ptr_rf => interp_eta1_rf
  interp_eta2_ptr_rf => interp_eta2_rf

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

!!$   ! jacobian times distribution function is stored
!!$       do j=1,num_pts2
!!$          do i=1, num_pts1
!!$             y = m%x2_cell(i,j)
!!$             x = m%x1_cell(i,j)
!!$             jac = m%jacobians_c(i,j)
!!$             data_out(i,j) = &
!!$                  jac / (2*sll_pi*init_obj%sigma_x*init_obj%sigma_y)*exp(-0.5_f64*( &
!!$                  (x-init_obj%xc)**2/init_obj%sigma_x**2 + &
!!$                  (y-init_obj%yc)**2/init_obj%sigma_y**2))
!!$          end do
!!$       end do

  print*, 'write mesh and distribution function'
 
  call write_scalar_field_2d(df)!,multiply_by_jacobian=.true.) 

  ! Initialize Linrood plan
  call new_linrood_plan(linrood, df)

  Print*, 'checking advection of a Gaussian in a uniform field'

  ! define uniform field on coarse_mesh (using stream function)
  call initialize_scalar_field_2d( &
       uniform_field, &
       "uniform_field", &
       m, &
       CELL_CENTERED_FIELD, &
       interp_eta1_ptr_sf, &
       interp_eta2_ptr_sf )

  ! components of field
  alpha1 = 1.0_f64
  alpha2 = 1.0_f64
  do i1 = 1, nc_eta1 
     do i2 = 1, nc_eta2
        FIELD_2D_AT_I( uniform_field, i1, i2 ) = alpha1 * m%x2_cell(i1,i2) &
             - alpha2 * m%x1_cell(i1,i2)
     end do
  end do
 
  print*, 'checking advection in rotating field' 
  ! define rotating field
  call initialize_scalar_field_2d( &
       rotating_field, &
       "rotating_field", &
       m, &
       CELL_CENTERED_FIELD, &
       interp_eta1_ptr_rf, &
       interp_eta2_ptr_rf )

  do i1 = 1, nc_eta1
     do i2 = 1, nc_eta2
        FIELD_2D_AT_I( rotating_field, i1, i2 ) = 0.5_f64*(m%x1_cell(i1,i2)**2 &
             + m%x2_cell(i1,i2)**2)
     end do
  end do

 print*, 'checking advection in incompressible swirling deformation field' 
  ! define incompressible field
  call initialize_scalar_field_2d(incompressible_field, "incompressible_field", &
       m, CELL_CENTERED_FIELD, interp_eta1_ptr_rf, interp_eta2_ptr_rf )
  do i1 = 1, nc_eta1
     do i2 = 1, nc_eta2
        FIELD_2D_AT_I( incompressible_field, i1, i2 ) = &
            (cos(sll_pi*2.0_f64*m%x1_cell(i1,i2)) &
             *cos(sll_pi*2.0_f64*m%x2_cell(i1,i2)) &
             -cos(sll_pi*2.0_f64*m%x1_cell(i1,i2)) & 
             -cos(sll_pi*2.0_f64*m%x2_cell(i1,i2))) & 
             /sll_pi/4.0_f64
     end do
  end do

  deltat = m%delta_eta1*.1_f64

  nbiter = 600
  do istep = 1, nbiter
     call linrood_step(linrood, df, incompressible_field, 0._8, deltat)
 
     call write_scalar_field_2d(df)!,multiply_by_jacobian=.true.)
     call write_scalar_field_2d(incompressible_field) 
  end do
  print *, 'Successful, exiting program.'

end program unit_test
