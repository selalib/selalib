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
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_timer
  implicit none

  character(len=256) :: filename
  type(sll_simulation_4d_vp_eulerian_cartesian_finite_volume)      :: simulation
  type(sll_logical_mesh_2d), pointer      :: mx,mv
  type(sll_time_mark)  :: t0 
  class(sll_coordinate_transformation_2d_base),pointer      :: tx,tv
  sll_real64, dimension(1:11) :: landau_params
  sll_real64 :: time

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

#define NCELL1 16
#define NCELL2 16
#define NCELL3 16
#define NCELL4 16
!!$#define NCELL1 16
!!$#define NCELL2 16
!!$#define NCELL3 32
!!$#define NCELL4 32
!!$!transport
!!$#define ETA1MIN -1.0_f64
!!$#define ETA1MAX 1.0_f64
!!$#define ETA2MIN -1.0_f64
!!$#define ETA2MAX 1.0_f64
!!$#define ETA3MIN -1.0_f64
!!$#define ETA3MAX 1.0_f64
!!$#define ETA4MIN -1.0_f64
!!$#define ETA4MAX 1.0_f64
!landau 1d sur xvx or 2 streams
!!$#define ETA1MIN -6.0_f64
!!$#define ETA1MAX 6.0_f64
!!$#define ETA2MIN -0.5_f64
!!$#define ETA2MAX 0.5_f64
!!$#define ETA3MIN 0.0_f64
!!$#define ETA3MAX 2.0_f64*sll_pi/0.2
!!$#define ETA4MIN 0.0_f64
!!$#define ETA4MAX 1.0_f64
!galaxy 1d sur xvx 
!!$#define ETA1MIN -4.0_f64
!!$#define ETA1MAX 4.0_f64
!!$#define ETA2MIN -4.0_f64
!!$#define ETA2MAX 4.0_f64
!!$#define ETA3MIN -12.0_f64
!!$#define ETA3MAX 12.0_f64
!!$#define ETA4MIN -12.0_f64
!!$#define ETA4MAX 12.0_f64
!!$!landau 1d sur yvy
!!$#define ETA1MIN -0.5_f64
!!$#define ETA1MAX 0.5_f64
!!$#define ETA2MIN -6.0_f64
!!$#define ETA2MAX 6.0_f64
!!$#define ETA3MIN 0.0_f64
!!$#define ETA3MAX 1.0_f64
!!$#define ETA4MIN 0.0_f64
!!$#define ETA4MAX 4.0_f64*sll_pi
!!$!landau 2D
#define ETA1MIN -6.0_f64
#define ETA1MAX 6.0_f64
#define ETA2MIN -6.0_f64
#define ETA2MAX 6.0_f64
#define ETA3MIN 0.0_f64
#define ETA3MAX 4.0_f64*sll_pi
#define ETA4MIN 0.0_f64
#define ETA4MAX 4.0_f64*sll_pi


#define TINI 0.0_f64
#define TMAX 6.1_f64
!#define TMAX 0._f64
#define CFL 0.7_f64
#define ELECMAX 1.0_f64 ! upper bound estimate for the electric field
#define EPSILON 0.005
#define TEST 5
! 0: x transport 1: landau damping 1d xvx  2: vx-transport
! 3: vy transport 4: y transport 5: landau 2d
!6: transport x-vx 7: transport y-vy 8: transport 2d
!9: landau damping 1d sur y-vy
!10: two-streams instability
!11: galaxy 1D test case
!12: galaxy 2D test case

#define DEG  2 ! polynomial degree
#define SCHEME 1
!0 Euler 1: Rung-Kutta 2 order 2:Rung-Kutta 4 order


  ! logical mesh for space coordinates
  mx => new_logical_mesh_2d( NCELL3, NCELL4 , &
       eta1_min=ETA3MIN, eta1_max=ETA3MAX, &
       eta2_min=ETA4MIN,eta2_max=ETA4MAX )


  ! logical mesh for velocity coordinates
  mv => new_logical_mesh_2d( NCELL1, NCELL2, &
       eta1_min=ETA1MIN, eta1_max=ETA1MAX, &
       eta2_min=ETA2MIN,eta2_max=ETA2MAX )


  tx => new_coordinate_transformation_2d_analytic( &
    'mapxy',          &
    mx,        &
    identity_x1,        &
    identity_x2,        &
    identity_jac11,       &
    identity_jac12,       &
    identity_jac21,       &
    identity_jac22, &
    (/0.0_f64,0.0_f64/) )

  tv => new_coordinate_transformation_2d_analytic( &
    'mapvxvy',          &
    mv,        &
    identity_x1,        &
    identity_x2,        &
    identity_jac11,       &
    identity_jac12,       &
    identity_jac21,       &
    identity_jac22, &
    (/0.0_f64,0.0_f64/) )


  ! define the values of the parameters for the landau initializer
    landau_params(1)=ETA3MIN
    landau_params(2)=ETA3MAX
    landau_params(3)=ETA4MIN
    landau_params(4)=ETA4MAX
    landau_params(5)= EPSILON  !epsilon in the landau
    landau_params(6)= DEG  ! polynomial interpolation degree
    landau_params(7)=CFL
    landau_params(8)=TEST
    landau_params(9)=ELECMAX
    landau_params(10)=SCHEME
    landau_params(11)=TINI
  ! initialize simulation object with the above parameters
    if(TEST==0) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_test_x_transport_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==1) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_landau_1d_xvx_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==2) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_test_vx_transport_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==3) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_test_vy_transport_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==4) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_test_y_transport_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==5) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_landau_2d_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==6) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_test_xvx_transport_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==7) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_test_yvy_transport_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==9) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_landau_1d_yvy_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==10) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_twostream_1d_xvx_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
    else if (TEST==11) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_galaxy_1d_xvx_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )
  else if (TEST==12) then
       call initialize_vp4d( &
            simulation, &
            mx,mv,tx,tv, &
            sll_galaxy_2d_initializer_v1v2x1x2, &
            landau_params, &
            TMAX )

    end if
  print *, 'Start time mark t0'
  call sll_set_time_mark(t0)
  call simulation%run( )
  time = sll_time_elapsed_since(t0)
  print *, 'time of simulation est  : ',time
  call sll_delete(simulation)

  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()


end program vp_cartesian_4d
