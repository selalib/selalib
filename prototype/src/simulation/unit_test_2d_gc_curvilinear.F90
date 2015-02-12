program guiding_center_2d_curvilinear
 use sll_simulation_2d_guiding_center_curvilinear_module
 implicit none
class(sll_simulation_base_class), pointer :: sim
character(len=256) :: filename
character(len=256) :: filename_local

  call get_command_argument(1, filename)
  if (len_trim(filename) == 0)then
    sim => new_guiding_center_2d_curvilinear( )
  else
    filename_local = trim(filename)
    sim => new_guiding_center_2d_curvilinear( filename_local )
  endif
print *,'#PASSED NEW GUIDING CENTER 2D CURVILINEAR'
call sim%run( )
print *,'#PASSED RUN'

end program guiding_center_2d_curvilinear