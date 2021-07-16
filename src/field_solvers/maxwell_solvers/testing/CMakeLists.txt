IF(NOT PGI)
  ADD_EXECUTABLE(test_maxwell_1d_fem test_maxwell_1d_fem.F90)
  TARGET_LINK_LIBRARIES(test_maxwell_1d_fem sll_maxwell_solvers)
  
  ADD_TEST(NAME maxwell_1d_fem COMMAND test_maxwell_1d_fem)
  SET_TESTS_PROPERTIES(maxwell_1d_fem PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  
  ADD_EXECUTABLE(test_maxwell_1d_fem_strongA test_maxwell_1d_fem_strongA.F90)
  TARGET_LINK_LIBRARIES(test_maxwell_1d_fem_strongA sll_maxwell_solvers)
  
  ADD_TEST(NAME maxwell_1d_fem_strongA COMMAND test_maxwell_1d_fem_strongA)
  SET_TESTS_PROPERTIES(maxwell_1d_fem_strongA PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  
  ADD_EXECUTABLE(test_maxwell_1d_ps test_maxwell_1d_ps.F90)
  TARGET_LINK_LIBRARIES(test_maxwell_1d_ps sll_maxwell_solvers)
  
  ADD_TEST(NAME maxwell_1d_ps COMMAND test_maxwell_1d_ps)
  SET_TESTS_PROPERTIES(maxwell_1d_ps PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ENDIF(NOT PGI)

ADD_EXECUTABLE(test_maxwell_2d_fdtd test_maxwell_2d_fdtd.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_fdtd sll_maxwell_solvers)

ADD_TEST(NAME maxwell_2d_fdtd COMMAND test_maxwell_2d_fdtd)
SET_TESTS_PROPERTIES(maxwell_2d_fdtd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_LIBRARY(maxwell_helper_functions STATIC m_maxwell_helper_functions.F90)

TARGET_LINK_LIBRARIES(maxwell_helper_functions sll_constants)
ADD_EXECUTABLE(test_maxwell_2d_diga_periodic test_maxwell_2d_diga_periodic.F90)
ADD_TEST(NAME maxwell_2d_diga_periodic COMMAND test_maxwell_2d_diga_periodic)
SET_TESTS_PROPERTIES(maxwell_2d_diga_periodic 
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_periodic maxwell_helper_functions sll_maxwell_solvers)
ADD_EXECUTABLE(test_maxwell_2d_diga_discrete test_maxwell_2d_discrete.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_discrete maxwell_helper_functions sll_maxwell_solvers)

ADD_EXECUTABLE(test_maxwell_2d_diga_wave test_maxwell_2d_diga_wave.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_diga_wave maxwell_helper_functions sll_maxwell_solvers)
ADD_TEST(NAME maxwell_2d_diga_wave COMMAND test_maxwell_2d_diga_wave)
SET_TESTS_PROPERTIES(maxwell_2d_diga_wave
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_maxwell_2d_pstd test_maxwell_2d_pstd.F90 )
TARGET_LINK_LIBRARIES(test_maxwell_2d_pstd sll_maxwell_solvers)
ADD_TEST(NAME maxwell_2d_pstd COMMAND test_maxwell_2d_pstd)
SET_TESTS_PROPERTIES(maxwell_2d_pstd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_maxwell_3d_pstd test_maxwell_3d_pstd.F90 )
TARGET_LINK_LIBRARIES(test_maxwell_3d_pstd sll_maxwell_solvers)
ADD_TEST(NAME maxwell_3d_pstd COMMAND test_maxwell_3d_pstd)
SET_TESTS_PROPERTIES(maxwell_3d_pstd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")



ADD_EXECUTABLE(test_maxwell_3d_fem_fft test_maxwell_3d_fem_fft.F90)
TARGET_LINK_LIBRARIES(test_maxwell_3d_fem_fft sll_maxwell_solvers)

ADD_TEST(NAME maxwell_3d_fem_fft COMMAND test_maxwell_3d_fem_fft)
SET_TESTS_PROPERTIES(maxwell_3d_fem_fft PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


ADD_EXECUTABLE(test_maxwell_2d_fem_fft test_maxwell_2d_fem_fft.F90)
TARGET_LINK_LIBRARIES(test_maxwell_2d_fem_fft sll_maxwell_solvers)

ADD_TEST(NAME maxwell_2d_fem_fft COMMAND test_maxwell_2d_fem_fft)
SET_TESTS_PROPERTIES(maxwell_2d_fem_fft PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


ADD_EXECUTABLE(test_maxwell_1d_fem_sm test_maxwell_1d_fem_sm.F90)
TARGET_LINK_LIBRARIES(test_maxwell_1d_fem_sm sll_maxwell_solvers)

ADD_TEST(NAME maxwell_1d_fem_sm COMMAND test_maxwell_1d_fem_sm)
SET_TESTS_PROPERTIES(maxwell_1d_fem_sm PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_maxwell_3d_fem test_maxwell_3d_fem.F90)
TARGET_LINK_LIBRARIES(test_maxwell_3d_fem sll_maxwell_solvers)

ADD_TEST(NAME maxwell_3d_fem COMMAND test_maxwell_3d_fem)
SET_TESTS_PROPERTIES(maxwell_3d_fem PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 120)


 ADD_EXECUTABLE(test_maxwell_1d_trafo test_maxwell_1d_trafo.F90)
 TARGET_LINK_LIBRARIES(test_maxwell_1d_trafo sll_maxwell_solvers )
 
 ADD_TEST(NAME maxwell_1d_trafo COMMAND test_maxwell_1d_trafo)
 SET_TESTS_PROPERTIES(maxwell_1d_trafo PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
 
 ADD_EXECUTABLE(test_maxwell_3d_trafo test_maxwell_3d_trafo.F90)
 TARGET_LINK_LIBRARIES(test_maxwell_3d_trafo sll_maxwell_solvers )
 
 ADD_TEST(NAME maxwell_3d_trafo COMMAND test_maxwell_3d_trafo)
 SET_TESTS_PROPERTIES(maxwell_3d_trafo PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 120) 

 ADD_EXECUTABLE(test_maxwell_clamped_1d_fem_sm test_maxwell_clamped_1d_fem_sm.F90)
 TARGET_LINK_LIBRARIES(test_maxwell_clamped_1d_fem_sm sll_maxwell_solvers)
 
 ADD_TEST(NAME maxwell_clamped_1d_fem_sm COMMAND test_maxwell_clamped_1d_fem_sm)
 SET_TESTS_PROPERTIES(maxwell_clamped_1d_fem_sm PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
 
  ADD_EXECUTABLE(test_maxwell_clamped_3d_fem test_maxwell_clamped_3d_fem.F90)
 TARGET_LINK_LIBRARIES(test_maxwell_clamped_3d_fem sll_maxwell_solvers)
 
 ADD_TEST(NAME maxwell_clamped_3d_fem COMMAND test_maxwell_clamped_3d_fem)
 SET_TESTS_PROPERTIES(maxwell_clamped_3d_fem PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 120)


 ADD_EXECUTABLE(test_maxwell_clamped_1d_trafo test_maxwell_clamped_1d_trafo.F90)
 TARGET_LINK_LIBRARIES(test_maxwell_clamped_1d_trafo sll_maxwell_solvers)
 
 ADD_TEST(NAME maxwell_clamped_1d_trafo COMMAND test_maxwell_clamped_1d_trafo)
 SET_TESTS_PROPERTIES(maxwell_clamped_1d_trafo PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
 
 ADD_EXECUTABLE(test_maxwell_clamped_3d_trafo test_maxwell_clamped_3d_trafo.F90)
 TARGET_LINK_LIBRARIES(test_maxwell_clamped_3d_trafo sll_maxwell_solvers )

 ADD_TEST(NAME maxwell_clamped_3d_trafo COMMAND test_maxwell_clamped_3d_trafo)
 SET_TESTS_PROPERTIES(maxwell_clamped_3d_trafo PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 120)
