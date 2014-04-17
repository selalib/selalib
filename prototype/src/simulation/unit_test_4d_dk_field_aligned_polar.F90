! Sample computation with the following characteristics:
! - drift-kinetic
! - 4D : r,\theta,z and v
! - parallel
! - filed aligned

program drift_kinetic_field_aligned_polar
  use sll_simulation_4d_drift_kinetic_field_aligned_polar_module
  use sll_collective
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_drift_kinetic_field_aligned_polar) :: simulation
  call sll_boot_collective()
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#Booting parallel environment...'
  endif

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call getarg(1, filename)
  filename_local = trim(filename)
  call simulation%init_from_file(filename_local)
  call simulation%run( )
  call delete(simulation)
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#reached end of dk4d_polar test'
    print *, '#PASSED'
  endif
  call sll_halt_collective()


end program drift_kinetic_field_aligned_polar
