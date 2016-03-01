# CTest
add_executable( test_pic_particles test_pic_particles.F90 )
target_link_libraries( test_pic_particles sll_pic_opt2d_particles )
add_test( NAME pic_particles COMMAND test_pic_particles )
set_tests_properties( pic_particles PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
