program guiding_center_2d_curvilinear
 use sll_simulation_2d_guiding_center_curvilinear_module
 implicit none
class(sll_simulation_base_class), pointer :: sim

print *,'#PASSED'
sim => new_guiding_center_2d_curvilinear()
print *,'#PASSED NEW GUIDING CENTER 2D CURVILINEAR'
call sim%run( )

print *,'#PASSED RUN'

end program guiding_center_2d_curvilinear