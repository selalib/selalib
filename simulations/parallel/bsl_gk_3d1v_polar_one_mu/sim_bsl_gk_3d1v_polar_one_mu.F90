! Sample computation with the following characteristics:
! - drift-kinetic with one mu
! -for the moment copy of drift_kinetic_polar
! - 4D : r,\theta,z and v
! - parallel

program sim_bsl_gk_3d1v_polar_one_mu
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   use sll_m_collective, only: &
      sll_s_boot_collective, &
      sll_f_get_collective_rank, &
      sll_s_halt_collective, &
      sll_v_world_collective

   use sll_m_sim_bsl_gk_3d1v_polar_one_mu, only: &
      sll_s_delete_dk4d_polar, &
      sll_t_simulation_4d_drift_kinetic_polar_one_mu

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   character(len=256) :: filename
   character(len=256) :: filename_local
   type(sll_t_simulation_4d_drift_kinetic_polar_one_mu) :: simulation
   call sll_s_boot_collective()
   if (sll_f_get_collective_rank(sll_v_world_collective) == 0) then
      print *, '#Booting parallel environment...'
   end if

   ! In this test, the name of the file to open is provided as a command line
   ! argument.
   call get_command_argument(1, filename)
   filename_local = trim(filename)
   call simulation%init_from_file(filename_local)
   call simulation%run()
   call sll_s_delete_dk4d_polar(simulation)
   if (sll_f_get_collective_rank(sll_v_world_collective) == 0) then
      print *, '#reached end of dk4d_polar test'
      print *, '#PASSED'
   end if
   call sll_s_halt_collective()

end program sim_bsl_gk_3d1v_polar_one_mu
