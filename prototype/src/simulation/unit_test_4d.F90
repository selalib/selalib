! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D cartesian: x, y, vx, vy (or x1, x2, x3, x4)
! - parallel

program vlasov_poisson_4d
  use sll_simulation_4d_vlasov_poisson_cartesian
  use sll_collective
  implicit none

  type(sll_simulation_4d_vlasov_poisson_cart) :: simulation

  print *, 'Booting parallel environment...'
  call sll_boot_collective() ! Wrap this up in something else
  ! call simulation%read_from_file("whatever") or something similar here
  call simulation%run( )

  print *, 'reached end of vp4d test'
  print *, 'PASSED'
  call delete_vp4d_par_cart(simulation)

  call sll_halt_collective()

end program vlasov_poisson_4d


