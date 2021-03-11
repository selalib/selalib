if(NOT PGI)
  add_executable(test_maxwell_1d_fem test_maxwell_1d_fem.F90)
  target_link_libraries(test_maxwell_1d_fem sll_maxwell_solvers)

  add_test(NAME maxwell_1d_fem COMMAND test_maxwell_1d_fem)
  set_tests_properties(maxwell_1d_fem PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED")
endif(NOT PGI)

add_executable(test_maxwell_2d_fdtd test_maxwell_2d_fdtd.F90)
target_link_libraries(test_maxwell_2d_fdtd sll_maxwell_solvers)

add_test(NAME maxwell_2d_fdtd COMMAND test_maxwell_2d_fdtd)
set_tests_properties(maxwell_2d_fdtd PROPERTIES PASS_REGULAR_EXPRESSION
                                                "PASSED")

add_library(maxwell_helper_functions STATIC m_maxwell_helper_functions.F90)

target_link_libraries(maxwell_helper_functions sll_constants)
add_executable(test_maxwell_2d_diga_periodic test_maxwell_2d_diga_periodic.F90)
add_test(NAME maxwell_2d_diga_periodic COMMAND test_maxwell_2d_diga_periodic)
set_tests_properties(maxwell_2d_diga_periodic PROPERTIES PASS_REGULAR_EXPRESSION
                                                         "PASSED")

target_link_libraries(test_maxwell_2d_diga_periodic maxwell_helper_functions
                      sll_maxwell_solvers)
add_executable(test_maxwell_2d_diga_discrete test_maxwell_2d_discrete.F90)
target_link_libraries(test_maxwell_2d_diga_discrete maxwell_helper_functions
                      sll_maxwell_solvers)

add_executable(test_maxwell_2d_diga_wave test_maxwell_2d_diga_wave.F90)
target_link_libraries(test_maxwell_2d_diga_wave maxwell_helper_functions
                      sll_maxwell_solvers)
add_test(NAME maxwell_2d_diga_wave COMMAND test_maxwell_2d_diga_wave)
set_tests_properties(maxwell_2d_diga_wave PROPERTIES PASS_REGULAR_EXPRESSION
                                                     "PASSED")

add_executable(test_maxwell_2d_pstd test_maxwell_2d_pstd.F90)
target_link_libraries(test_maxwell_2d_pstd sll_maxwell_solvers)
add_test(NAME maxwell_2d_pstd COMMAND test_maxwell_2d_pstd)
set_tests_properties(maxwell_2d_pstd PROPERTIES PASS_REGULAR_EXPRESSION
                                                "PASSED")

add_executable(test_maxwell_3d_pstd test_maxwell_3d_pstd.F90)
target_link_libraries(test_maxwell_3d_pstd sll_maxwell_solvers)
add_test(NAME maxwell_3d_pstd COMMAND test_maxwell_3d_pstd)
set_tests_properties(maxwell_3d_pstd PROPERTIES PASS_REGULAR_EXPRESSION
                                                "PASSED")
