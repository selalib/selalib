program guiding_center_2d_cartesian
 use sll_simulation_2d_guiding_center_cartesian_module
 implicit none
class(sll_simulation_base_class), pointer :: sim

sim => new_guiding_center_2d_cartesian()
call sim%run( )



sim => new_guiding_center_2d_cartesian()
call sim%run( )


 print *,'#PASSED'

end program guiding_center_2d_cartesian
