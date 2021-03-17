add_executable(test_lobalap_analytic test_lobalap_analytic.F90)
target_link_libraries(test_lobalap_analytic
                      sll_general_coordinate_elliptic_solvers)

add_executable(test_lobalap_discrete test_lobalap_discrete.F90)
target_link_libraries(test_lobalap_discrete
                      sll_general_coordinate_elliptic_solvers)
