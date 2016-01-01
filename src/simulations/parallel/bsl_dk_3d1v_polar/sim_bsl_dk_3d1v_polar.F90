! Sample computation with the following characteristics:
! - drift-kinetic
! - 4D : r,\theta,z and v
! - parallel

program sim_bsl_dk_3d1v_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_sim_bsl_dk_3d1v_polar, only: &
    sll_s_delete_dk4d_polar, &
    sll_t_simulation_4d_drift_kinetic_polar

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_since, &
    sll_t_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_t_simulation_4d_drift_kinetic_polar) :: simulation
  type(sll_t_time_mark)  :: t0
  sll_real64 :: time 

  call sll_s_boot_collective()
  if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
    print *, '#Booting parallel environment...'
    print *, '#Start time mark t0'
    call sll_s_set_time_mark(t0)
  endif

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call get_command_argument(1, filename)
  filename_local = trim(filename)
  call simulation%init_from_file(filename_local)
  call simulation%run( )
  call sll_s_delete_dk4d_polar(simulation)
  if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
    print *, '#reached end of dk4d_polar test'
    time = sll_f_time_elapsed_since(t0)
    print *, '#time elapsed since t0 : ',time
    print *, '#PASSED'
  endif
  call sll_s_halt_collective()


end program sim_bsl_dk_3d1v_polar
