add_executable(test_particle_sort test_particle_sort.F90)
target_link_libraries(test_particle_sort sll_random_deviate_generators
                      sll_pic_opt_utilities sll_pic_opt_particles sll_meshes)
add_test(NAME pic_opt_particle_sort COMMAND test_particle_sort)
set_tests_properties(pic_opt_particle_sort PROPERTIES PASS_REGULAR_EXPRESSION
                                                      "PASSED")

add_executable(test_pic_utilities test_pic_utilities.F90)
target_link_libraries(test_pic_utilities sll_pic_opt_utilities sll_meshes)
add_test(NAME pic_opt_utilities COMMAND test_pic_utilities)
set_tests_properties(pic_opt_utilities PROPERTIES PASS_REGULAR_EXPRESSION
                                                  "PASSED")

add_executable(test_charge_to_density test_charge_to_density.F90)
target_link_libraries(test_charge_to_density sll_pic_opt_utilities sll_meshes)
add_test(NAME pic_opt_charge_to_density COMMAND test_charge_to_density)
set_tests_properties(pic_opt_charge_to_density
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
