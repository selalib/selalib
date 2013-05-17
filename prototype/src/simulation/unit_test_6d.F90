! Sample computation with the following characteristics:
! - vlasov-poisson
! - 6D cartesian: x, y, z, vx, vy, vz (or x1, x2, x3, x4, x5, x6)
! - parallel

program vlasov_poisson_6d
  use sll_simulation_6d_vlasov_poisson_cartesian
  use sll_collective
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_6d_vlasov_poisson_cart) :: simulation

  print *, 'Booting parallel environment for 6D simulation...'
  call sll_boot_collective() ! Wrap this up in something else

  print *, 'Proceed to run simulation.'
  call flush()

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call getarg(1, filename)
  filename_local = trim(filename)
  call simulation%init_from_file(filename_local)
  call simulation%run( )
  call delete_vp6d_par_cart(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'
  call sll_halt_collective()

end program vlasov_poisson_6d


