add_executable(test_pic_poisson_2d test_pic_poisson_2d.F90)
target_link_libraries(test_pic_poisson_2d sll_pic_field_solvers
                      sll_particle_mesh_coupling sll_poisson_solvers)

add_test(NAME pic_poisson_2d COMMAND test_pic_poisson_2d)
set_tests_properties(pic_poisson_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
