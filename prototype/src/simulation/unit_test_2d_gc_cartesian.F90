program guiding_center_2d_cartesian
  use sll_simulation_2d_guiding_center_cartesian_module
  implicit none
  class(sll_simulation_base_class), pointer :: sim
  character(len=256) :: filename
  character(len=256) :: filename_local

  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => new_guiding_center_2d_cartesian( )
  else
    filename_local = trim(filename)
    sim => new_guiding_center_2d_cartesian( filename_local )
  endif
  call sim%run( )
  print *,'#PASSED'

end program guiding_center_2d_cartesian
