! Sample computation with the following characteristics:
! - drift kinetic
! - 4D: x1,x2,x3,v3 (or v1) with cartesian coordinate 
! - parallel

program dk_cartesian_4d
#include "sll_working_precision.h"
  use sll_simulation_4d_drift_kinetic_cartesian_finite_volume
  use sll_collective
  use sll_constants
  use sll_logical_meshes
  use sll_common_array_initializers_module
  implicit none

  type(sll_simulation_4d_drift_kinetic_cart_finite_volume)      :: simulation
  type(sll_logical_mesh_4d), pointer      :: mx
  sll_real64, dimension(1:2) :: landau_params

  print *, 'Booting parallel environment...'
  call sll_boot_collective() ! Wrap this up somewhere else

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  !call getarg(1, filename)
  !filename_local = trim(filename)

  ! To initialize the simulation type, there should be two options. One is to
  ! initialize from a file:
  
  !call simulation%init_from_file(filename_local)

  
  ! The second is to initialize 'manually' with a routine whose parameters
  ! allow to configure the different types of objects in the simulation. For
  ! instance, the type of coordinate mapping. Here we use both methods while
  ! we develop and sort out the interfaces.
  ! Eventually, when using the module, one should only need to use one 
  ! way to initialize the simulation object, in development we are using them
  ! both...

! hardwired, this should be consistent with whatever is read from a file
#define NCELL1 32
#define NCELL2 32
#define NCELL3 4
#define NCELL4 8

  ! logical mesh for space coordinates
  mx => new_logical_mesh_4d( NCELL1, NCELL2,NCELL3, NCELL4 , &
       eta1_min=-6.0_f64, eta1_max=6.0_f64,eta2_min=0.0_f64,eta2_max=1.0_f64)




  ! define the values of the parameters for the landau initializer
  landau_params(1) = 0.1
  landau_params(2) = 2.0*sll_pi

  ! initialize simulation object with the above parameters
  call initialize_dk4d( &
       simulation, &
       mx, &
       sll_landau_initializer_dk_test_4d, &
       landau_params )

  call simulation%run( )
  call sll_delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()


end program dk_cartesian_4d


