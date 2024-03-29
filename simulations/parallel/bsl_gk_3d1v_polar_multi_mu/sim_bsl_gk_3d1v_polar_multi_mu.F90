! Sample computation with the following characteristics:
! - drift-kinetic with several mu
! comes from 4d_drift_kinetic_polar_one_mu
! - 4D : r,\theta,z and v
! - parallel

program sim_bsl_gk_3d1v_polar_multi_mu
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_collective, only: &
      sll_s_boot_collective, &
      sll_f_get_collective_rank, &
      sll_s_halt_collective, &
      sll_v_world_collective

   use sll_m_sim_bsl_gk_3d1v_polar_multi_mu, only: &
      sll_s_delete_dk4d_polar, &
      sll_t_simulation_4d_drift_kinetic_polar_multi_mu

   use sll_m_timer, only: &
      sll_s_set_time_mark, &
      sll_f_time_elapsed_since, &
      sll_t_time_mark

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   character(len=256) :: filename
   character(len=256) :: filename_local
   type(sll_t_simulation_4d_drift_kinetic_polar_multi_mu) :: simulation
   type(sll_t_time_mark)  :: t0
   sll_real64 :: time!,time1,time2
   call sll_s_boot_collective()
   if (sll_f_get_collective_rank(sll_v_world_collective) == 0) then
      print *, '#Booting parallel environment...'
      print *, '#Start time mark t0'
      call sll_s_set_time_mark(t0)
   end if

   ! In this test, the name of the file to open is provided as a command line
   ! argument.
   call get_command_argument(1, filename)
   filename_local = trim(filename)
   call simulation%init_from_file(filename_local)
!if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
!  time1 = sll_f_time_elapsed_since(t0)
!  print *, '#time elapsed for init : ',time1
!endif
   call simulation%run()
!if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
!  time2 = sll_f_time_elapsed_since(t0)
!  print *, '#time elapsed for run : ', time2 - time1
!endif
   call sll_s_delete_dk4d_polar(simulation)
   if (sll_f_get_collective_rank(sll_v_world_collective) == 0) then
      print *, '#reached end of dk4d_polar_multi_mu test'
      time = sll_f_time_elapsed_since(t0)
      print *, '#time elapsed since t0 : ', time
      print *, '#PASSED'
   end if
   call sll_s_halt_collective()

end program sim_bsl_gk_3d1v_polar_multi_mu
