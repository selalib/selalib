program unit_test

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
  use numeric_constants
  use distribution_function
  use sll_module_mapped_meshes_2d
  use user_geometry_functions
  use geometry_functions
  use sll_scalar_field_initializers_base
  use sll_landau_2d_initializer
  use sll_module_interpolators_1d_base
  use sll_cubic_spline_interpolator_1d
  use sll_advection_field
  implicit none

  sll_int32 :: nc_eta1, nc_eta2
  type(sll_mapped_mesh_2d_analytic), target           :: mesh2d
  class(sll_mapped_mesh_2d_base), pointer             :: m
  type(sll_distribution_function_2d)                  :: df 
  character(32)                                       :: name = 'dist_func'
  character(len=4)                                    :: cstep
  type(init_landau_2d), target                        :: init_landau
  class(scalar_field_2d_initializer_base), pointer    :: p_init_f
  type(cubic_spline_1d_interpolator)                  :: spline_interp
  type(hamiltonian_advection_field_2d)                :: hamiltonian
  sll_int32  :: ierr, istep
  sll_int32  :: ix, iv, nnode_x1, nnode_v1
  sll_real64, dimension(:), allocatable               :: feet_x1
  sll_real64, dimension(:), allocatable               :: feet_x2
  sll_real64, dimension(:), allocatable               :: field_x1
  sll_real64, dimension(:), allocatable               :: field_x2
  nc_eta1 = 64
  nc_eta2 = 64

  SLL_ALLOCATE( feet_x1(nc_eta1+1), ierr )
  SLL_ALLOCATE( feet_x2(nc_eta2+1), ierr )
  SLL_ALLOCATE( field_x1(nc_eta1+1), ierr )
  SLL_ALLOCATE( field_x2(nc_eta2+1), ierr )

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
  m => mesh2d
  print*, 'initialization of distribution_function'
  call init_landau%initialize(0.001_f64)
  p_init_f => init_landau
  print *, 'completed initialization of distribution function'

  call initialize_distribution_function_2d( &
       df, &
       1.0_f64, &
       1.0_f64, &
       name, &
       m, &
       NODE_CENTERED_FIELD, &
       p_init_f )
  print*, 'write mesh and distribution function'

  istep = 0
  call int2string(istep,cstep)
  df%name = trim(name)//cstep
  call write_scalar_field_2d(df,multiply_by_jacobian=.true.) 

  print *, 'Initialize cubic spline interpolator'
  call spline_interp%initialize( nc_eta1+1, 0.0_f64, 1.0_f64, PERIODIC_SPLINE )

  print *, 'Initialize the hamiltonian advection field'
  call initialize_advection_field_2d( hamiltonian, 1.0_f64, 1.0_f64, &
       "hamiltonian", mesh2d, NODE_CENTERED_FIELD )

  print *, 'Do an advection step'
  ! Advection step consists of:
  ! 1. Build the spline interpolant for the distribution function line
  ! 2. For each node in the grid,
  call spline_interp%compute_interpolants( df%data(:,iv) )
  
#if 0
  do iv=1,nc_eta2+1
       do ix=1,nc_eta1+1
        ! is this always done in the eta1,eta2 space? The value of eta1_min
        ! is hardwired as 0.0 here. Should this change for generality?
        feet_x(ix) 
#endif

  print *, 'Successful, exiting program.'

end program unit_test
