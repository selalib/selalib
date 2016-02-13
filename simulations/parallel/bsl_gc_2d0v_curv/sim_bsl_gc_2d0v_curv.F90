program sim_bsl_gc_2d0v_curv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_v_world_collective

  use sll_m_sim_bsl_gc_2d0v_curv, only: &
    sll_s_delete_guiding_center_2d_curvilinear, &
    sll_f_new_guiding_center_2d_curvilinear, &
    sll_t_simulation_2d_guiding_center_curvilinear

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_since, &
    sll_t_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !class(sll_c_simulation_base_class), pointer :: sim
  class(sll_t_simulation_2d_guiding_center_curvilinear), pointer :: sim  
  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_t_time_mark)  :: t0
  sll_real64 :: time
  sll_int32 :: count
  sll_int32 :: i
  sll_int32 :: num_min
  sll_int32 :: num_max
  character(len=256) :: str






  call sll_s_boot_collective()
  if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
    print *, '#Start time mark t0'
    call sll_s_set_time_mark(t0)
    print *, '#Booting parallel environment...'
  endif

  count = command_argument_count()
  call get_command_argument(1, filename)

  if (len_trim(filename) == 0)then
    sim => sll_f_new_guiding_center_2d_curvilinear( )
  else
    filename_local = trim(filename)
    call get_command_argument(2, str)
    if(len_trim(str) == 0)then
      sim => sll_f_new_guiding_center_2d_curvilinear( filename_local )
      call sim%run( )
    else
      read(str , *) num_max
      num_min = 0
      call get_command_argument(3, str)
      if(len_trim(str) .ne. 0)then
        num_min = num_max
        read(str , *) num_max
      endif
      do i=num_min,num_max
        sim => sll_f_new_guiding_center_2d_curvilinear( filename_local, i)
        call sim%run( )
        call sll_s_delete_guiding_center_2d_curvilinear( sim )
      enddo  
    endif
  endif

  if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
    print *, '#reached end of guiding_center_2d_curvilinear test'
    time = sll_f_time_elapsed_since(t0)
    print *, '#time elapsed since t0 : ',time
    print *, '#PASSED'
  endif


 



end program sim_bsl_gc_2d0v_curv
