program sim_bsl_gc_2d0v_cart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_sim_bsl_gc_2d0v_cart, only: &
    sll_f_new_guiding_center_2d_cartesian

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  class(sll_c_simulation_base_class), pointer :: sim
  character(len=256) :: filename
  character(len=256) :: filename_local

  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => sll_f_new_guiding_center_2d_cartesian( )
  else
    filename_local = trim(filename)
    sim => sll_f_new_guiding_center_2d_cartesian( filename_local )
  endif
  call sim%run( )
  print *,'#PASSED'

end program sim_bsl_gc_2d0v_cart
