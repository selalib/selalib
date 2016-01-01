program sim_bsl_ad_2d0v_curv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_sim_bsl_ad_2d0v_curv, only: &
    sll_s_delete_analytic_field_2d_curvilinear, &
    sll_f_new_analytic_field_2d_curvilinear, &
    sll_t_simulation_2d_analytic_field_curvilinear

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  class(sll_t_simulation_2d_analytic_field_curvilinear), pointer :: sim  
  character(len=256) :: filename
  character(len=256) :: filename_local
  sll_int32 :: count
  sll_int32 :: i
  sll_int32 :: num_min
  sll_int32 :: num_max
  character(len=256) :: str

  count = command_argument_count()
 
  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => sll_f_new_analytic_field_2d_curvilinear( )
  else
    filename_local = trim(filename)
    call get_command_argument(2, str)
    if(len_trim(str) == 0)then
      sim => sll_f_new_analytic_field_2d_curvilinear( filename_local )
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
        sim => sll_f_new_analytic_field_2d_curvilinear( filename_local, i)
        call sim%run( )
        call sll_s_delete_analytic_field_2d_curvilinear( sim )
      enddo  
    endif  
  endif
  
  print *,'#PASSED'




end program sim_bsl_ad_2d0v_curv
