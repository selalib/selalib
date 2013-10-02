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
!VG!  use sll_common_array_initializers_module
  implicit none

  sll_int32 :: world_size
  sll_int32 :: my_rank

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_DK_hybrid)  :: simulation
  type(sll_logical_mesh_4d), pointer :: logical_mesh4d
  

  ! Parallelization initialization
  print *, 'Booting parallel environment...'
  call sll_boot_collective() ! Wrap this up somewhere else
  world_size = sll_get_collective_size(sll_world_collective)
  my_rank    = sll_get_collective_rank(sll_world_collective)

  ! Reading of the input file 'sim4d_DK_hybrid_input.txt'
  call getarg(1,filename)
  filename_local = trim(filename)
  call simulation%init_from_file(filename_local)

  ! logical mesh for space coordinates
  logical_mesh4D => new_logical_mesh_4d( &
    simulation%nc_x1,simulation%nc_x2, &
    simulation%nc_x3,simulation%nc_x4,eta1_min=0.1_f64, &
    eta1_max=14.5_f64,eta2_min=0.0_f64,eta2_max=2.0*sll_pi, &
    eta3_min=0._f64,eta3_max=1508._f64,eta4_min=-6._f64,eta4_max=6._f64)

  ! initialize 4D drift-kinetic Vlasov
  call initialize(simulation, &
    world_size, &
    my_rank, &
    logical_mesh4D)

  !** Fisrt step ***
  call first_step_4d_DK_hybrid(simulation)

!VG!  call simulation%run( )
!VG!  call delete(simulation)

  print *, 'reached end of 4d DK hybrid test'
  print *, 'PASSED'

  call sll_halt_collective()
end program DK_hybrid_4d


