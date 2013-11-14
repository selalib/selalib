! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D cartesian: x,y, vx, vy (or x1, x2)
! - parallel for vlasov
! - sequential for poisson 

program vlasov_poisson_4d
  use sll_simulation_4d_vlasov_parallel_poisson_sequential_cartesian
  use sll_collective
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  class(sll_simulation_base_class), pointer :: sim
  call sll_boot_collective()
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#Booting parallel environment...'
  endif

  
  ! In this test, the name of the file to open is provided as a command line
  ! argument.

  
  call getarg(1, filename)
  filename_local = trim(filename)
  
  sim => new_vlasov_par_poisson_seq_cart( filename_local )
  
  call sim%run()
  
  !call simulation%init_from_file(filename_local)
  !call simulation%run( )
  !call delete_vp2d_par_cart(simulation)
  
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#reached end of sim4d_vp_cart test'
    print *, '#PASSED'
  endif
  call sll_halt_collective()


end program vlasov_poisson_4d
