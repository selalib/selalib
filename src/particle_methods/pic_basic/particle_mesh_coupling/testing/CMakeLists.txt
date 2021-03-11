add_executable(test_particle_mesh_coupling_spline_1d
               test_particle_mesh_coupling_spline_1d.F90)
target_link_libraries(test_particle_mesh_coupling_spline_1d
                      sll_particle_mesh_coupling sll_particle_groups)
add_test(NAME particle_mesh_coupling_spline_1d
         COMMAND test_particle_mesh_coupling_spline_1d)
set_tests_properties(particle_mesh_coupling_spline_1d
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_particle_mesh_coupling_spline_2d
               test_particle_mesh_coupling_spline_2d.F90)
target_link_libraries(test_particle_mesh_coupling_spline_2d
                      sll_particle_mesh_coupling sll_particle_groups)
add_test(NAME particle_mesh_coupling_spline_2d
         COMMAND test_particle_mesh_coupling_spline_2d)
set_tests_properties(particle_mesh_coupling_spline_2d
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
