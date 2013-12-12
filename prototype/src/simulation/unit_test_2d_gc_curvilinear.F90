program guiding_center_2d_curvilinear
 use sll_simulation_2d_guiding_center_curvilinear_module
 implicit none
class(sll_simulation_base_class), pointer :: sim


sim => new_guiding_center_2d_curvilinear()
call sim%run( )


 print *,'#PASSED'

end program guiding_center_2d_curvilinear