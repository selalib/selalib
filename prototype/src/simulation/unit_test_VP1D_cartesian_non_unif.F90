! Sample computation with the following characteristics:
! - vlasov-poisson
! - 1Dx1D cartesian: x, vx (or x1, x2)
! - non uniform mesh
! - sequential


program unit_test_VP1D_cartesian_non_unif
  use simulation_VP1D_cartesian_non_unif
  implicit none

  type(sll_simulation_VP1D_cartesian_non_unif) :: simulation

  print *, '#Begin of VP1D_cartesian_non_unif test'
  call simulation%run( )

  print *, '#reached end of VP1D_cartesian_non_unif test'
  print *, '#PASSED'
  !call delete_vp4d_par_cart(simulation)


end program unit_test_VP1D_cartesian_non_unif


