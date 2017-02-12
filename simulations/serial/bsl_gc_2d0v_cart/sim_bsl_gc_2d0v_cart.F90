program sim_bsl_gc_2d0v_cart
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_sim_bsl_gc_2d0v_cart, only: &
     initialize_guiding_center_2d_cartesian, &
     sll_simulation_2d_guiding_center_cartesian, &
    sll_f_new_guiding_center_2d_cartesian

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  type(sll_simulation_2d_guiding_center_cartesian) :: sim
  character(len=256) :: filename
  character(len=256) :: filename_local

  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
     call initialize_guiding_center_2d_cartesian(sim)
  else
    filename_local = trim(filename)
     call initialize_guiding_center_2d_cartesian(sim, filename_local )
  endif
  call sim%run( )
  print *,'#PASSED'

end program sim_bsl_gc_2d0v_cart
