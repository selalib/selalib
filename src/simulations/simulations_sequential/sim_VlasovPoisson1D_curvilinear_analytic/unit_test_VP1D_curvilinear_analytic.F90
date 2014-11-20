! Sample computation with the following characteristics:
! - vlasov-poisson
! - 1Dx1D cartesian: x, vx (or x1, x2)
! - curvilinear mesh with analytic mapping
! - sequential


program unit_test_VP1D_curvilinear_analytic
  use simulation_VP1D_curvilinear_analytic
  implicit none

  type(sll_simulation_VP1D_curvilinear_analytic) :: simulation

  print *, '#Begin of VP1D_curvilinear_analytic test'
  call simulation%run( )

  print *, '#reached end of VP1D_curvilinear_analytic test'
  print *, '#PASSED'
  !call delete_vp4d_par_cart(simulation)


end program unit_test_VP1D_curvilinear_analytic


