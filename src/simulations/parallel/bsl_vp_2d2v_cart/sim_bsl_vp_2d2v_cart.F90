! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D cartesian: x, y, vx, vy (or x1, x2, x3, x4)
! - parallel

program sim_bsl_vp_2d2v_cart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_collective, only: &
    sll_boot_collective, &
    sll_halt_collective

  use sll_m_sim_bsl_vp_2d2v_cart, only: &
    delete_vp4d_par_cart, &
    sll_simulation_4d_vlasov_poisson_cart

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_vlasov_poisson_cart) :: simulation

  print *, 'Booting parallel environment...'
  call sll_boot_collective() ! Wrap this up somewhere else

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call get_command_argument(1, filename)
  filename_local = trim(filename)
  call simulation%init_from_file(filename_local)
  call simulation%run( )
  call delete_vp4d_par_cart(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()


end program sim_bsl_vp_2d2v_cart



