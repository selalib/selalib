program guiding_center_2d_curvilinear_mudpack
 use sll_simulation_2d_guiding_center_curvilinear_mudpack_module
 implicit none
 class(sll_simulation_2d_guiding_center_curvilinear_mudpack), pointer :: sim
 !type(sll_logical_mesh_2d), pointer      :: M
 !class(sll_coordinate_transformation_2d_base), pointer :: transformation
 
 

sim => new_guiding_center_2d_curvilinear_mudpack()
call sim%run( )


 print *,'#PASSED'

end program guiding_center_2d_curvilinear_mudpack
