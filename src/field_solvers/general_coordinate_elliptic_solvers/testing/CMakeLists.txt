add_executable(
  test_general_coordinates_elliptic_solver
  m_init_functions.F90 test_general_coordinates_elliptic_solver.F90)

target_link_libraries(test_general_coordinates_elliptic_solver
                      sll_general_coordinate_elliptic_solvers)

add_test(NAME general_coordinate_elliptic_solver
         COMMAND test_general_coordinates_elliptic_solver)

set_tests_properties(general_coordinate_elliptic_solver
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 2000)

add_executable(test_lobalap_analytic test_lobalap_analytic.F90)
target_link_libraries(test_lobalap_analytic
                      sll_general_coordinate_elliptic_solvers)


add_executable(test_lobalap_discrete test_lobalap_discrete.F90)
target_link_libraries(test_lobalap_discrete
                      sll_general_coordinate_elliptic_solvers)

add_executable(test_gces_full_periodic_prototype
               test_full_periodic_prototype.F90)
target_link_libraries(test_gces_full_periodic_prototype
                      sll_general_coordinate_elliptic_solvers)
