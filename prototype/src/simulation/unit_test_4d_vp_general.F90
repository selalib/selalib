! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with arbitrary coordinate transformation
!   in the x,y variables.
! - parallel

program vlasov_poisson_4d_general
  use sll_simulation_4d_vlasov_poisson_general
  use sll_collective
  use sll_constants
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_vp_general)      :: simulation
  type(sll_logical_mesh_2d), pointer      :: mx
  type(sll_logical_mesh_2d), pointer      :: mv
  class(sll_coordinate_transformation_2d_base), pointer :: transformation_x
  sll_real64, dimension(1:3) :: landau_params

  print *, 'Booting parallel environment...'
  call sll_boot_collective() ! Wrap this up somewhere else

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call getarg(1, filename)
  filename_local = trim(filename)

  ! To initialize the simulation type, there should be two options. One is to
  ! initialize from a file:
  
!  call simulation%init_from_file(filename_local)
  
  ! The second is to initialize 'manually' with a routine whose parameters
  ! allow to configure the different types of objects in the simulation. For
  ! instance, the type of coordinate mapping. Here we use both methods while
  ! we develop and sort out the interfaces.
  ! Eventually, when using the module, one should only need to use one 
  ! way to initialize the simulation object.

! hardwired, this should be consistent with whatever is read from a file
#define NPTS1 32
#define NPTS2 32

  ! logical mesh for space coordinates
  mx => sll_new_logical_mesh_2d( NPTS1, NPTS2 )

  ! logical mesh for velocity coordinates
  mv => sll_new_logical_mesh_2d( NPTS1, NPTS2, &
       eta1_min=-6.0_f64, eta1_max=6.0_f64, &
       eta2_min=-6.0_f64, eta2_max=6.0_f64)

  ! coordinate transformation associated with space coordinates
  transformation_x => new_coordinate_transformation_2d_analytic( &
       "analytic_identity_transformation", &
       mx, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22 )

  ! define the values of the parameters for the landau initializer
  landau_params(1) = 0.1
  landau_params(2) = 2.0*sll_pi
  landau_params(3) = 2.0*sll_pi

  ! initialize simulation object with the above parameters
  call initialize( &
       simulation, &
       mx, &
       mv, &
       transformation_x, &
       sll_landau_initializer_4d, &
       landau_params )

  call simulation%run( )
  call delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()


end program vlasov_poisson_4d_general


