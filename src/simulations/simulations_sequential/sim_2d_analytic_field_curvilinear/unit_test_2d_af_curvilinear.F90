program analytic_field_2d_curvilinear
  use sll_simulation_2d_analytic_field_curvilinear_module
  implicit none
  class(sll_simulation_base_class), pointer :: sim
  character(len=256) :: filename
  character(len=256) :: filename_local

  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => new_analytic_field_2d_curvilinear( )
  else
    filename_local = trim(filename)
    sim => new_analytic_field_2d_curvilinear( filename_local )
  endif
  call sim%run( )
  print *,'#PASSED'

end program analytic_field_2d_curvilinear
