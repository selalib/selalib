! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with arbitrary coordinate 
!   transformation
!   in the x,y variables.
! - parallel

program vlasov_poisson_4d_polar
#include "sll_working_precision.h"
#include "sll_constants.h"

  use sll_simulation_4d_vlasov_poisson_polar
  use sll_collective
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  implicit none

  character(len=256)                  :: filename
  character(len=256)                  :: filename_local
  type(sll_simulation_4d_vp_polar)    :: simulation
  type(sll_logical_mesh_2d), pointer  :: mx
  type(sll_logical_mesh_2d), pointer  :: mv
  sll_real64                          :: params(6)

  class(sll_coordinate_transformation_2d_base), pointer :: transformation

  call sll_boot_collective() ! Wrap this up somewhere else

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call getarg(1, filename)
  filename_local = trim(filename)

  ! To initialize the simulation type, there should be two options. One is to
  ! initialize from a file:
  
  call simulation%init_from_file(filename_local)
  
  ! The second is to initialize 'manually' with a routine whose parameters
  ! allow to configure the different types of objects in the simulation. For
  ! instance, the type of coordinate mapping. Here we use both methods while
  ! we develop and sort out the interfaces.
  ! Eventually, when using the module, one should only need to use one 
  ! way to initialize the simulation object, in development we are using them
  ! both...

! hardwired, this should be consistent with whatever is read from a file
#define NPTS1 32
#define NPTS2 32
#define NPTS3 32
#define NPTS4 32

  ! logical mesh for space coordinates
  mx => new_logical_mesh_2d( NPTS1, NPTS2,       & 
       eta1_min=.2_f64, eta1_max=.8_f64,         &
       eta2_min=.0_f64, eta2_max=2.0_f64*sll_pi)

  ! logical mesh for velocity coordinates
  mv => new_logical_mesh_2d( NPTS3, NPTS4, &
       eta1_min=-6.0_f64, eta1_max=6.0_f64, &
       eta2_min=-6.0_f64, eta2_max=6.0_f64)

  ! coordinate transformation associated with space coordinates
  transformation => new_coordinate_transformation_2d_analytic( &
       "analytic_polar_transformation", &
       mx, &
       polar_x1, &
       polar_x2, &
       polar_jac11, &
       polar_jac12, &
       polar_jac21, &
       polar_jac22 )

  ! initialize simulation object with the above parameters
  call initialize_vp4d_polar( &
       simulation, &
       mx, &
       mv, &
       transformation, &
       sll_periodic_gaussian_initializer_4d, &
       params )

!  function defined in  parallel_array_initializers/sll_common_array_initializers_module.F90
!  sll_periodic_gaussian_initializer_4d(x,y,vx,xy) = 
!  val = alpha*exp(-0.5_f64*((x -xc )**2+(y -yc )**2)) + &
!        beta *exp(-0.5_f64*((vx-vxc)**2+(vy-vyc)**2))

  params(1) = 0.5 !xc
  params(2) = 0.0 !yc
  params(3) = 0.0 !vxc
  params(4) = 0.0 !vyc
  params(5) = 1.0 !alpha
  params(6) = 0.0 !beta


  call simulation%run( )
  call delete(simulation)

  print *, 'PASSED'

  call sll_halt_collective()


end program vlasov_poisson_4d_polar
