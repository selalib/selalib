add_executable( test_particle_group_2d test_particle_group_2d.F90 )
target_link_libraries( test_particle_group_2d sll_pic_opt_particles)
add_test( NAME pic_opt_particle_group_2d COMMAND test_particle_group_2d )
set_tests_properties( pic_opt_particle_group_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

add_executable( test_particle_group_4d test_particle_group_4d.F90 )
target_link_libraries( test_particle_group_4d sll_pic_opt_particles )
add_test( NAME pic_opt_particle_group_4d COMMAND test_particle_group_4d )
set_tests_properties( pic_opt_particle_group_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
