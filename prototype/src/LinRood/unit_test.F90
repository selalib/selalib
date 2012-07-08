program unit_test
#include "sll_working_precision.h"
#include "sll_field_2d.h"
#include "sll_memory.h"
  use numeric_constants
  use geometry_functions
  use sll_module_mapped_meshes_2d
  use sll_gaussian_2d_initializer
  use distribution_function
  use sll_scalar_field_2d
  use sll_linrood
  implicit none

  sll_int32 :: nc_eta1, nc_eta2
  type(sll_mapped_mesh_2d_analytic), target :: mesh2d
  class(sll_mapped_mesh_2d_base), pointer   :: m
  type(sll_distribution_function_2d)   :: df 
  character(32)  :: name = 'dist_func'
  character(len=4) :: cstep
  type(init_gaussian_2d), target :: init_gaussian
  class(scalar_field_2d_initializer_base), pointer    :: p_init_f
  type(scalar_field_2d) :: uniform_field, rotating_field 
  type(linrood_plan) :: linrood
  
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
  m => mesh2d

  print*, 'initialization of distribution_function'

  call init_gaussian%initialize( CELL_CENTERED_FIELD, 0.5_f64, 0.5_f64, 0.1_f64, 0.1_f64 )
  p_init_f => init_gaussian

  call initialize_distribution_function_2d( &
       df, &
       1.0_f64, &
       1.0_f64, &
       name, &
       m, &
       CELL_CENTERED_FIELD, &
       p_init_f )
  print*, 'write mesh and distribution function'

  istep = 0
  call int2string(istep,cstep)
  df%name = trim(name)//cstep
  
  call write_scalar_field_2d(df)!,multiply_by_jacobian=.true.) 

  ! Initialize Linrood plan
  call new_linrood_plan(linrood, df)

  Print*, 'checking advection of a Gaussian in a uniform field'

  ! define uniform field on coarse_mesh (using stream function)
  call initialize_scalar_field_2d(uniform_field, "uniform_field", m, CELL_CENTERED_FIELD)
  ! components of field
  alpha1 = 0.0_f64
  alpha2 = 10.0_f64
  do i1 = 1, nc_eta1 
     do i2 = 1, nc_eta2
        FIELD_2D_AT_I( uniform_field, i1, i2 ) = alpha1 * m%x2_cell(i1,i2) &
             - alpha2 * m%x1_cell(i1,i2)
     end do
  end do
 
  print*, 'checking advection in rotating field' 
  ! define rotating field
  call initialize_scalar_field_2d(rotating_field, "rotating_field", m, CELL_CENTERED_FIELD)
  do i1 = 1, nc_eta1
     do i2 = 1, nc_eta2
        FIELD_2D_AT_I( rotating_field, i1, i2 ) = 0.5_f64*(m%x1_cell(i1,i2)**2 &
             + m%x2_cell(i1,i2)**2)
     end do
  end do


  print *, 'Successful, exiting program.'

end program unit_test
