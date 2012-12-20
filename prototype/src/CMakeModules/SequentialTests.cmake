ADD_TEST(NAME memory COMMAND test_memory)
SET_TESTS_PROPERTIES(memory PROPERTIES TIMEOUT 20)

IF(CMAKE_BUILD_TYPE STREQUAL DEBUG)
   ADD_TEST(NAME assert COMMAND test_assert)
   SET(passRegex "Assertion error triggered in file")
   SET_TESTS_PROPERTIES(assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}")
ENDIF(CMAKE_BUILD_TYPE STREQUAL DEBUG)

ADD_TEST(NAME constants COMMAND test_constants)


ADD_TEST(NAME timer COMMAND test_timer)

ADD_TEST(NAME tridiagonal COMMAND test_tridiagonal)
ADD_TEST(NAME lagrange COMMAND test_lagrange)

ADD_TEST(NAME toeplitz_penta_diagonal COMMAND test_toeplitz_penta_diagonal)
SET_TESTS_PROPERTIES(toeplitz_penta_diagonal PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME newton_raphson COMMAND test_newton_raphson)

ADD_TEST(NAME splines COMMAND test_splines ) 
SET_TESTS_PROPERTIES(splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME splines_arbitrary_degree COMMAND test_arbitrary_degree_splines)
SET_TESTS_PROPERTIES(splines_arbitrary_degree PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME quintic_splines COMMAND test_quintic_splines)
SET_TESTS_PROPERTIES(quintic_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME odd_degree_splines COMMAND test_odd_degree_splines)
SET_TESTS_PROPERTIES(odd_degree_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME cubic_non_uniform_splines COMMAND test_non_unif_splines)
SET_TESTS_PROPERTIES(cubic_non_uniform_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


ADD_TEST(NAME integration COMMAND test_integration)

IF(FFTPACK_ENABLED)
  ADD_TEST(NAME periodic_interp COMMAND test_periodic_interp)
ENDIF()

ADD_TEST(NAME fft COMMAND test_fft)

IF(NOT STDF95)
   ADD_TEST(NAME utilities COMMAND test_utilities)
   IF(FFTPACK_ENABLED)
      ADD_TEST(NAME poisson_solvers COMMAND test_poisson_1d)
   ENDIF(FFTPACK_ENABLED)

   ADD_TEST(NAME poisson_3d_periodic_seq COMMAND test_poisson_3d_periodic_seq)
   SET_TESTS_PROPERTIES(poisson_3d_periodic_seq 
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ADD_TEST(NAME qns2d_with_finite_diff_seq 
            COMMAND test_qns2d_with_finite_diff_seq)
   SET_TESTS_PROPERTIES(qns2d_with_finite_diff_seq 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ADD_TEST(NAME qns2d_angular_spectral_method_seq 
            COMMAND test_qns2d_angular_spectral_method_seq)
   SET_TESTS_PROPERTIES(qns2d_angular_spectral_method_seq 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ADD_TEST(NAME interpolators COMMAND test_interpolators_1d test_interpolators_2d)

   ADD_TEST(NAME fields COMMAND test_scalar_field)
   ADD_TEST(NAME time_splitting COMMAND test_time_splitting)
   ADD_TEST(NAME distribution_function COMMAND test_distribution_function)
   ADD_TEST(NAME advection_field COMMAND test_advection_field)

   IF(FFTW_FOUND)
      ADD_TEST(NAME maxwell_2d_pstd 
            COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_maxwell_2d_pstd)
      SET_TESTS_PROPERTIES(maxwell_2d_pstd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ENDIF()

ENDIF()

ADD_TEST(NAME WENO COMMAND test_WENO_interp test_WENO_recon)
ADD_TEST(NAME quintic_interpolators_1d COMMAND test_quintic_interpolators_1d)
SET_TESTS_PROPERTIES(quintic_interpolators_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME quintic_interpolators_1d_nonuniform COMMAND test_quintic_interpolators_1d_nonuniform)
SET_TESTS_PROPERTIES(quintic_interpolators_1d_nonuniform PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME odd_degree_interpolators_1d COMMAND test_odd_degree_interpolators_1d)
SET_TESTS_PROPERTIES(odd_degree_interpolators_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME odd_degree_interpolators_1d_nonuniform COMMAND test_odd_degree_interpolators_1d_nonuniform)
SET_TESTS_PROPERTIES(odd_degree_interpolators_1d_nonuniform PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME electric_field_accumulators COMMAND test_e_field_accumulator_2d)

ADD_TEST(NAME mapped_meshes COMMAND test_mapped_meshes_1d
				    test_mapped_meshes_2d)

ADD_TEST(NAME ode_solvers COMMAND test_implicit_ode_nonuniform)

ADD_TEST(NAME BSL COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_bsl_2d)
SET_TESTS_PROPERTIES(BSL PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME low_level_file_io COMMAND test_io)
SET_TESTS_PROPERTIES(low_level_file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME maxwell_2d_fdtd 
         COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_maxwell_2d_fdtd)
SET_TESTS_PROPERTIES(maxwell_2d_fdtd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

IF(FORTRANCL_FOUND)
   ADD_TEST(NAME opencl COMMAND test_opencl)
   SET_TESTS_PROPERTIES(opencl PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF(FORTRANCL_FOUND)

IF(MUDPACK_ENABLED AND Fortran_COMPILER STREQUAL "GFORTRAN")
   ADD_TEST(NAME mudpack COMMAND tmud34sp tmud24sp )
   SET(mudpack PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()
