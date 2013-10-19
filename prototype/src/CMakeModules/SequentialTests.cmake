ADD_TEST(NAME memory COMMAND test_memory)
SET_TESTS_PROPERTIES(memory PROPERTIES TIMEOUT 20)

IF(CMAKE_BUILD_TYPE STREQUAL DEBUG)
   ADD_TEST(NAME assert COMMAND test_assert)
   SET(passRegex "Assertion error triggered in file")
   SET_TESTS_PROPERTIES(assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}")
ENDIF(CMAKE_BUILD_TYPE STREQUAL DEBUG)

ADD_TEST(NAME constants                 COMMAND test_constants)
ADD_TEST(NAME timer                     COMMAND test_timer)
SET_TESTS_PROPERTIES(timer PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME logical_meshes            COMMAND test_logical_meshes)
ADD_TEST(NAME tridiagonal               COMMAND test_tridiagonal)
ADD_TEST(NAME lagrange                  COMMAND test_lagrange)
ADD_TEST(NAME toeplitz_penta_diagonal   COMMAND test_toeplitz_penta_diagonal)
ADD_TEST(NAME newton_raphson            COMMAND test_newton_raphson)
ADD_TEST(NAME splines                   COMMAND test_splines) 
ADD_TEST(NAME splines_arbitrary_degree  COMMAND test_arbitrary_degree_splines)
ADD_TEST(NAME quintic_splines           COMMAND test_quintic_splines)
ADD_TEST(NAME odd_degree_splines        COMMAND test_odd_degree_splines)
ADD_TEST(NAME cubic_non_uniform_splines COMMAND test_non_unif_splines)
ADD_TEST(NAME integration               COMMAND test_integration)
ADD_TEST(NAME lagrange_interpolation    COMMAND test_lagrange_interpolation)

SET_TESTS_PROPERTIES(logical_meshes PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(toeplitz_penta_diagonal PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(splines_arbitrary_degree PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(quintic_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(odd_degree_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(cubic_non_uniform_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(integration PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(lagrange_interpolation PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME periodic_interp COMMAND test_periodic_interp)

ADD_TEST(NAME fft COMMAND test_fft)

IF(NOT STDF95)
   ADD_TEST(NAME utilities COMMAND test_utilities)
   ADD_TEST(NAME poisson_solvers COMMAND test_poisson_1d)

   ADD_TEST(NAME poisson_3d_periodic_seq COMMAND test_poisson_3d_periodic_seq)
   SET_TESTS_PROPERTIES(poisson_3d_periodic_seq 
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ADD_TEST(NAME qns2d COMMAND test_qn_solver_2d)
   SET_TESTS_PROPERTIES(qns2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
#consider merging the following 2 tests
   ADD_TEST(NAME interpolators COMMAND test_interpolators_1d test_interpolators_2d)
   ADD_TEST(NAME arb_deg_spline_interpolator COMMAND test_arb_deg_spline_interpolators_2d)
   ADD_TEST(NAME arb_deg_spline_interpolator_1d COMMAND test_arb_deg_spline_interpolators_1d)
   ADD_TEST(NAME fields COMMAND test_scalar_field)
   ADD_TEST(NAME time_splitting COMMAND test_time_splitting)
   ADD_TEST(NAME distribution_function COMMAND test_distribution_function)
   ADD_TEST(NAME advection_field COMMAND test_advection_field)
   ADD_TEST(NAME coordinate_transformations COMMAND test_coordinate_transformations_2d)
   ADD_TEST(NAME general_coordinate_elliptic_solver COMMAND test_general_coordinates_elliptic_solver)
   ADD_TEST(NAME charac COMMAND test_characteristics)

   SET_TESTS_PROPERTIES(coordinate_transformations PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ADD_TEST(NAME fields_2d_alternative COMMAND test_scalar_field_alternative)
   SET_TESTS_PROPERTIES(fields_2d_alternative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET_TESTS_PROPERTIES(general_coordinate_elliptic_solver PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(charac PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET_TESTS_PROPERTIES(arb_deg_spline_interpolator PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(arb_deg_spline_interpolator_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   IF(FFTW_ENABLED)
      ADD_TEST(NAME maxwell_2d_pstd COMMAND test_maxwell_2d_pstd)
      SET_TESTS_PROPERTIES(maxwell_2d_pstd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ENDIF()

ENDIF()

ADD_TEST(NAME WENO COMMAND test_WENO_interp test_WENO_recon)
ADD_TEST(NAME quintic_1d COMMAND test_quintic_interpolators_1d)
SET_TESTS_PROPERTIES(quintic_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME quintic_1d_nonuniform COMMAND test_quintic_interpolators_1d_nonuniform)
SET_TESTS_PROPERTIES(quintic_1d_nonuniform PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME odd_degree_1d COMMAND test_odd_degree_interpolators_1d)
SET_TESTS_PROPERTIES(odd_degree_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME odd_degree_1d_nonuniform COMMAND test_odd_degree_interpolators_1d_nonuniform)
SET_TESTS_PROPERTIES(odd_degree_1d_nonuniform PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME electric_field_accumulators COMMAND test_e_field_accumulator_2d)

IF(NOT STDF95)
ADD_TEST(NAME mapped_meshes COMMAND test_mapped_meshes_1d
				    test_mapped_meshes_2d)

ADD_TEST(NAME ode_solvers COMMAND test_implicit_ode_nonuniform)

ADD_TEST(NAME BSL COMMAND bsl_1d_cubic_uniform_periodic
                          bsl_1d_cubic_nonuniform_periodic
                          bsl_1d_cubic_uniform_compact
                          bsl_1d_cubic_nonuniform_compact
                          bsl_1d_quintic_uniform_compact
                          bsl_1d_quintic_nonuniform_compact)

SET_TESTS_PROPERTIES(BSL PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()


IF(NOT STDF95)
ADD_TEST(NAME maxwell_2d_fdtd COMMAND test_maxwell_2d_fdtd)
SET_TESTS_PROPERTIES(maxwell_2d_fdtd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()

IF(FORTRANCL_FOUND)
   ADD_TEST(NAME opencl COMMAND test_opencl)
   SET_TESTS_PROPERTIES(opencl PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF(FORTRANCL_FOUND)
