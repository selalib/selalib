! Sample computation with the following characteristics:
! - vlasov-poisson
! - 6D cartesian: x, y, z, vx, vy, vz (or x1, x2, x3, x4, x5, x6)
! - parallel

program vlasov_poisson_6d
  use sll_simulation_6d_vlasov_poisson_cartesian
  use sll_collective
  implicit none

  type(sll_simulation_6d_vlasov_poisson_cart) :: simulation

  print *, 'Booting parallel environment for 6D simulation...'
  call sll_boot_collective() ! Wrap this up in something else
  ! call simulation%read_from_file("whatever") or something similar here
  print *, 'Proceed to run simulation.'
  call flush()
  call simulation%run( )

  print *, 'reached end of vp4d test'
  print *, 'PASSED'
  call delete_vp6d_par_cart(simulation)

  call sll_halt_collective()

end program vlasov_poisson_6d


