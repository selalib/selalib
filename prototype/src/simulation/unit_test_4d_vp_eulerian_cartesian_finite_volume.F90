! Sample computation with the following characteristics:
! - drift kinetic
! - 4D: x1,x2,x3,v3 (or v1) with cartesian coordinate 
! - parallel

program vp_cartesian_4d
#include "sll_working_precision.h"
  use sll_simulation_4d_vp_eulerian_cartesian_finite_volume_module
  use sll_collective
  use sll_constants
  use sll_logical_meshes
  use sll_common_array_initializers_module
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_vp_eulerian_cartesian_finite_volume)      :: simulation
  type(sll_logical_mesh_4d), pointer      :: mx
  sll_real64, dimension(1:6) :: landau_params

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
#define NCELL1 4
#define NCELL2 4
#define NCELL3 32
#define NCELL4 32
#define ETA1MIN -6.0_f64
#define ETA1MAX 6.0_f64
#define ETA2MIN -6.0_f64
#define ETA2MAX 6.0_f64
#define ETA3MIN 0.0_f64
#define ETA3MAX 1.0_f64
#define ETA4MIN 0.0_f64
#define ETA4MAX 1.0_f64

#define DEG 2   ! polynomial degree


  ! logical mesh for space coordinates
  mx => new_logical_mesh_4d( NCELL1, NCELL2,NCELL3, NCELL4 , &
       eta1_min=ETA1MIN, eta1_max=ETA1MAX, &
       eta2_min=ETA2MIN,eta2_max=ETA2MAX , &
       eta3_min=ETA3MIN, eta3_max=ETA3MAX, &
       eta4_min=ETA4MIN,eta4_max=ETA4MAX )




  ! define the values of the parameters for the landau initializer
    landau_params(1)=ETA3MIN
    landau_params(2)=ETA3MAX
    landau_params(3)=ETA4MIN
    landau_params(4)=ETA4MAX
    landau_params(5)= 0.1
    landau_params(6)= DEG  ! polynomial interpolation degree

  ! initialize simulation object with the above parameters
  call initialize_vp4d( &
       simulation, &
       mx, &
       sll_landau_initializer_v1v2x1x2, &
       landau_params, &
       1.0_f64 )


  call simulation%run( )
  call delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()


end program vp_cartesian_4d


