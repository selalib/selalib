! Sample computation with the following characteristics:
! - vlasov-poisson
! - 6D cartesian: x, y, z, vx, vy, vz (or x1, x2, x3, x4, x5, x6)
! - parallel

program sim_bsl_vp_3d3v_cart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use iso_fortran_env, only: &
    output_unit

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_s_halt_collective

  use sll_m_sim_bsl_vp_3d3v_cart, only: &
    sll_s_delete_vp6d_par_cart, &
    sll_t_simulation_6d_vlasov_poisson_cart

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_t_simulation_6d_vlasov_poisson_cart) :: simulation

  print *, 'Booting parallel environment for 6D simulation...'
  call sll_s_boot_collective() ! Wrap this up in something else

  print *, 'Proceed to run simulation.'
  flush( output_unit )

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call get_command_argument(1, filename)
  filename_local = trim(filename)
  call simulation%init_from_file(filename_local)
  call simulation%run( )
  call sll_s_delete_vp6d_par_cart(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'
  call sll_s_halt_collective()

end program sim_bsl_vp_3d3v_cart


