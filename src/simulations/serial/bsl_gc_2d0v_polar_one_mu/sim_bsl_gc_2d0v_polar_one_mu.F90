program sim_bsl_gc_2d0v_polar_one_mu
  use sll_m_sim_bsl_gc_2d0v_polar_one_mu
  implicit none
  class(sll_simulation_base_class), pointer :: sim
  character(len=256) :: filename
  character(len=256) :: filename_local
  
  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => new_guiding_center_2d_polar_one_mu( )
  else
    filename_local = trim(filename)
    sim => new_guiding_center_2d_polar_one_mu( filename_local )
  endif
  call sim%run( )
  print *,'#PASSED'

end program sim_bsl_gc_2d0v_polar_one_mu
