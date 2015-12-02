! Sample computation with the following characteristics:
! - drift-kinetic with several mu
! comes from 4d_drift_kinetic_polar_one_mu
! - 4D : r,\theta,z and v
! - parallel

program sim_bsl_gk_3d1v_polar_multi_mu
#include "sll_working_precision.h"
  use sll_m_sim_bsl_gk_3d1v_polar_multi_mu
  use sll_m_collective
  use sll_m_timer
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_drift_kinetic_polar_multi_mu) :: simulation
  type(sll_time_mark)  :: t0
  sll_real64 :: time,time1,time2 
  call sll_boot_collective()
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#Booting parallel environment...'
    print *, '#Start time mark t0'
    call sll_set_time_mark(t0)
  endif

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call get_command_argument(1, filename)
  filename_local = trim(filename)
  call simulation%init_from_file(filename_local)
!if(sll_get_collective_rank(sll_world_collective)==0)then
!  time1 = sll_time_elapsed_since(t0)
!  print *, '#time elapsed for init : ',time1
!endif     
  call simulation%run( )
!if(sll_get_collective_rank(sll_world_collective)==0)then
!  time2 = sll_time_elapsed_since(t0)
!  print *, '#time elapsed for run : ', time2 - time1
!endif  
  call delete_dk4d_polar(simulation)
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#reached end of dk4d_polar_multi_mu test'
    time = sll_time_elapsed_since(t0)
    print *, '#time elapsed since t0 : ',time
    print *, '#PASSED'
  endif
  call sll_halt_collective()


end program sim_bsl_gk_3d1v_polar_multi_mu