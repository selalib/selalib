add_executable(test_particle_visualization_interface
               test_particle_visualization_interface.F90)
target_link_libraries(test_particle_visualization_interface
                      sll_pic_visualization sll_pic_sampling sll_pic_weights)

set(NML_FILE ${CMAKE_CURRENT_SOURCE_DIR}/distribution_for_the_test.nml)
add_test(NAME particle_visualization_interface
         COMMAND test_particle_visualization_interface ${NML_FILE})
set_tests_properties(particle_visualization_interface
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
