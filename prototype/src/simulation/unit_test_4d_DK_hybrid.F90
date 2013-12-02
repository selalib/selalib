! Sample computation with the following characteristics:
! - drift kinetic
! - 4D: x1,x2,x3,v3 (or v1) with cartesian coordinate 
! - parallel

program DK_hybrid_4d
#include "sll_working_precision.h"
  use sll_simulation_4d_DK_hybrid_module
  use sll_collective
  use sll_constants
  use sll_logical_meshes
  use sll_common_coordinate_transformations
  use sll_module_coordinate_transformations_2d
  implicit none

  sll_int32 :: world_size
  sll_int32 :: my_rank

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_DK_hybrid)  :: simulation
  type(sll_logical_mesh_4d), pointer :: logical_mesh4d
  type(sll_logical_mesh_2d), pointer :: logical_mesh2d
  class(sll_coordinate_transformation_2d_base), pointer :: transf_xy

  ! Parallelization initialization
  print *, 'Booting parallel environment...'
  call sll_boot_collective() ! Wrap this up somewhere else
  print*, 'init from file'
  world_size = sll_get_collective_size(sll_world_collective)
  my_rank    = sll_get_collective_rank(sll_world_collective)
  ! Reading of the input file 'sim4d_DK_hybrid_input.txt'
  call getarg(1,filename)
  filename_local = trim(filename)
  call simulation%init_from_file(filename_local)

  print*, 'get logical mesh 4d'
  !*** logical mesh for space coordinates ***
  logical_mesh4D => new_logical_mesh_4d( &
    simulation%nc_x1,simulation%nc_x2, &
    simulation%nc_x3,simulation%nc_x4,eta1_min=0.0_f64, &
    eta1_max=1._f64,eta2_min=0.0_f64,eta2_max=1._f64, &
    eta3_min=simulation%phi_min,eta3_max=simulation%phi_max, &
    eta4_min=simulation%vpar_min,eta4_max=simulation%vpar_max)
  
  logical_mesh2d => new_logical_mesh_2d( &
    simulation%nc_x1,simulation%nc_x2)

  !*** coordinate transformation associated with space coordinates ***
  transf_xy => new_coordinate_transformation_2d_analytic( &
       "polar_transformation", &
       logical_mesh2d, &
       x1_polar_f, &
       x2_polar_f, &
       deriv_x1_polar_f_eta1, &
       deriv_x1_polar_f_eta2, &
       deriv_x2_polar_f_eta1, &
       deriv_x2_polar_f_eta2, &
       (/0.1_f64,1.0_f64/)) ! these were the default values for this map

  !*** initialize 4D drift-kinetic Vlasov ***
  call initialize(simulation, &
    world_size, &
    my_rank, &
    logical_mesh4D, &
    transf_xy)

  if (my_rank.eq.0) &
    call transf_xy%write_to_file()

  !*** Fisrt step ***
  if (simulation%my_rank.eq.0) &
    print*,': ---> First step'
  call first_step_4d_DK_hybrid(simulation)

  !*** Global loop ***
  if (simulation%my_rank.eq.0) &
    print*,': ---> Run'
  call simulation%run( )

  !*** Erase the memory used for the simulation ***
  if (simulation%my_rank.eq.0) &
    print*,': ---> Delete'
  call delete(simulation)

  print *, 'reached end of 4d DK hybrid test'
  print *, 'PASSED'

  call sll_halt_collective()
end program DK_hybrid_4d


