ADD_EXECUTABLE(test_maxwell_1d_fem test_maxwell_1d_fem.F90)
TARGET_LINK_LIBRARIES(test_maxwell_1d_fem sll_maxwell_solvers)

ADD_TEST(NAME maxwell_1d_fem COMMAND test_maxwell_1d_fem)
SET_TESTS_PROPERTIES(maxwell_1d_fem PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_maxwell_2d_fdtd test_maxwell_2d_fdtd.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_maxwell_solvers)

ADD_TEST(NAME maxwell_2d_fdtd COMMAND test_maxwell_2d_fdtd)
SET_TESTS_PROPERTIES(maxwell_2d_fdtd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_maxwell_2d_diga_periodic test_maxwell_2d_diga_periodic.F90 m_maxwell_helper_functions.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_periodic sll_maxwell_solvers)

ADD_EXECUTABLE(test_maxwell_2d_diga_wave test_maxwell_2d_diga_wave.F90 m_maxwell_helper_functions.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_wave sll_maxwell_solvers)

IF(HDF5_IS_PARALLEL)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)

ADD_EXECUTABLE(test_maxwell_2d_diga_discrete test_maxwell_2d_discrete.F90 m_maxwell_helper_functions.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_discrete sll_maxwell_solvers)

ADD_TEST(NAME maxwell_2d_diga_periodic COMMAND test_maxwell_2d_diga_periodic)
SET_TESTS_PROPERTIES(maxwell_2d_diga_periodic 
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME maxwell_2d_diga_wave COMMAND test_maxwell_2d_diga_wave)
SET_TESTS_PROPERTIES(maxwell_2d_diga_wave
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_maxwell_2d_pstd test_maxwell_2d_pstd.F90 )
TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_maxwell_solvers)
ADD_TEST(NAME maxwell_2d_pstd COMMAND test_maxwell_2d_pstd)
SET_TESTS_PROPERTIES(maxwell_2d_pstd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
