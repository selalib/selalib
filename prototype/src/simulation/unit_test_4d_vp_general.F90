! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D cartesian: x, y, vx, vy (or x1, x2, x3, x4)
! - parallel

program vlasov_poisson_4d_general
  use sll_simulation_4d_vlasov_poisson_general
  use sll_collective
  use sll_module_mapped_meshes_2d
  use geometry_functions
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_vp_general)      :: simulation
  class(sll_mapped_mesh_2d_base), pointer :: map2d

  print *, 'Booting parallel environment...'
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
  ! way to initialize the simulation object.

! hardwired, this should be consistent with whatever is read from a file
#define NPTS1 32
#define NPTS2 32
  map2d => new_mesh_2d_analytic( &
       "map_a", &
       NPTS1, &
       NPTS2, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22 )

  call initialize( simulation, map2d )

  call simulation%run( )
  call delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()


end program vlasov_poisson_4d_general


