program sim_bsl_dk_3d1v_curv_field_aligned
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_boot_collective, &
    sll_get_collective_rank, &
    sll_halt_collective, &
    sll_world_collective

  use sll_m_sim_bsl_dk_3d1v_curv_field_aligned, only: &
    sll_t_sim_sl_dk_3d1v_curv_field_aligned

  use sll_m_timer, only: &
    sll_set_time_mark, &
    sll_time_elapsed_since, &
    sll_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_sim_sl_dk_3d1v_curv_field_aligned) :: sim
  character(len=256)                            :: filename
  type(sll_time_mark)  :: t0
  sll_real64 :: time 

  call sll_boot_collective()
  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#Booting parallel environment...'
    print *, '#Start time mark t0'
    call sll_set_time_mark(t0)
  endif
  call get_command_argument(1, filename)

  if (len_trim(filename) == 0)then
    print *,'#please give input filename'
    stop
  else
    call sim%init_from_file(trim(filename))  
    call sim%run()
  endif

  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *, '#reached end of sim_sl_dk_3d1v_curv_field_aligned'
    time = sll_time_elapsed_since(t0)
    print *, '#time elapsed since t0 : ',time
    print *, '#PASSED'
  endif
  call sll_halt_collective()

  

end program sim_bsl_dk_3d1v_curv_field_aligned