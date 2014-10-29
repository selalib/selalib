ADD_TEST(NAME memory COMMAND test_memory)
SET_TESTS_PROPERTIES(memory PROPERTIES TIMEOUT 20)

IF(CMAKE_BUILD_TYPE STREQUAL DEBUG)
   ADD_TEST(NAME assert COMMAND test_assert)
   SET(passRegex "Assertion error triggered in file")
   SET_TESTS_PROPERTIES(assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}")
ENDIF(CMAKE_BUILD_TYPE STREQUAL DEBUG)

ADD_TEST(NAME constants                 COMMAND test_constants)
ADD_TEST(NAME timer                     COMMAND test_timer)
ADD_TEST(NAME tridiagonal               COMMAND test_tridiagonal)
ADD_TEST(NAME lagrange                  COMMAND test_lagrange)
ADD_TEST(NAME toeplitz_penta_diagonal   COMMAND test_toeplitz_penta_diagonal)
ADD_TEST(NAME cubic_splines             COMMAND test_splines) 
ADD_TEST(NAME splines_arbitrary_degree  COMMAND test_arbitrary_degree_splines)
ADD_TEST(NAME cubic_non_uniform_splines COMMAND test_non_unif_splines)
ADD_TEST(NAME integration               COMMAND test_integration)
ADD_TEST(NAME lagrange_interpolation    COMMAND test_lagrange_interpolation)
ADD_TEST(NAME hermite_interpolation     COMMAND test_hermite_interpolation)
ADD_TEST(NAME pic_particles             COMMAND test_pic_particles)
ADD_TEST(NAME pic_initializers          COMMAND test_pic_initializers)
ADD_TEST(NAME pic_accumulator           COMMAND test_pic_accumulator)
ADD_TEST(NAME pic_particle_sort         COMMAND test_particle_sort_2d)

ADD_TEST(NAME deboor_spline             COMMAND test_deboor)

SET_TESTS_PROPERTIES(timer PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(toeplitz_penta_diagonal PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(cubic_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(splines_arbitrary_degree PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(deboor_spline PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(cubic_non_uniform_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(integration PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(lagrange_interpolation PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(hermite_interpolation PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(pic_particles PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(pic_initializers PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(pic_accumulator PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(pic_particle_sort PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME fft COMMAND test_fft)

   ADD_TEST(NAME reduction COMMAND test_reduction)
   SET_TESTS_PROPERTIES(reduction PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ADD_TEST(NAME utilities COMMAND test_utilities)
   ADD_TEST(NAME poisson_solvers COMMAND test_poisson_1d)


#   ADD_TEST(NAME poisson_2d_dirichlet_cartesian COMMAND 
#     test_poisson_2d_dirichlet_cartesian )
#   SET_TESTS_PROPERTIES( poisson_2d_dirichlet_cartesian PROPERTIES 
#     PASS_REGULAR_EXPRESSION "PASSED" )


   ADD_TEST(NAME poisson_3d_periodic_seq COMMAND test_poisson_3d_periodic_seq)
   SET_TESTS_PROPERTIES(poisson_3d_periodic_seq 
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

#consider merging the following 2 tests
   ADD_TEST(NAME fields COMMAND test_scalar_field)
   ADD_TEST(NAME time_splitting COMMAND test_time_splitting)
   ADD_TEST(NAME distribution_function COMMAND test_distribution_function)
   ADD_TEST(NAME advection_field COMMAND test_advection_field)
   ADD_TEST(NAME coordinate_transformations COMMAND test_coordinate_transformations_2d)
   ADD_TEST(NAME fields_2d_alternative COMMAND test_scalar_field_alternative)
   ADD_TEST(NAME fields_1d_alternative COMMAND test_scalar_fields_1d_alternative)	

IF(PYTHON3_FOUND)

   ADD_TEST(NAME coordinate_transformation_multipatch_2d 
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMAND test_coordinate_transformation_multipatch_2d)
   SET_TESTS_PROPERTIES(coordinate_transformation_multipatch_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ADD_TEST(NAME scalar_field_multipatch_2d
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMAND test_scalar_field_multipatch_2d)
   SET_TESTS_PROPERTIES(scalar_field_multipatch_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


ELSE()

   MESSAGE(STATUS "python3 not found")

ENDIF(PYTHON3_FOUND)

   SET_TESTS_PROPERTIES(reduction PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ADD_TEST(NAME characteristics_1d_explicit_euler COMMAND test_characteristics_1d_explicit_euler)
   ADD_TEST(NAME characteristics_1d_explicit_euler_conservative
     COMMAND test_characteristics_1d_explicit_euler_conservative)
   ADD_TEST(NAME characteristics_1d_trapezoid COMMAND test_characteristics_1d_trapezoid)
   ADD_TEST(NAME characteristics_1d_trapezoid_conservative
     COMMAND test_characteristics_1d_trapezoid_conservative)
   ADD_TEST(NAME characteristics_2d_explicit_euler COMMAND test_characteristics_2d_explicit_euler)
   ADD_TEST(NAME characteristics_2d_explicit_euler_conservative 
     COMMAND test_characteristics_2d_explicit_euler_conservative)
   ADD_TEST(NAME characteristics_2d_verlet COMMAND test_characteristics_2d_verlet)
   ADD_TEST(NAME advection_1d_periodic COMMAND test_advection_1d_periodic)
   ADD_TEST(
     NAME
     advection_1d_non_uniform_cubic_splines 
     COMMAND
     test_advection_1d_non_uniform_cubic_splines)
   ADD_TEST(NAME advection_1d_BSL COMMAND test_advection_1d_BSL)
   ADD_TEST(NAME advection_1d_CSL COMMAND test_advection_1d_CSL)
   ADD_TEST(NAME advection_1d_PSM COMMAND test_advection_1d_PSM)
   ADD_TEST(NAME advection_2d_BSL COMMAND test_advection_2d_BSL)
   ADD_TEST(NAME advection_2d_CSL COMMAND test_advection_2d_CSL)
   ADD_TEST(NAME advection_2d_tensor_product COMMAND test_advection_2d_tensor_product)
   ADD_TEST(NAME gyroaverage_polar_hermite COMMAND test_gyroaverage_2d_polar_hermite)
   ADD_TEST(NAME gyroaverage_polar_splines COMMAND test_gyroaverage_2d_polar_splines)
   ADD_TEST(NAME gyroaverage_polar_pade COMMAND test_gyroaverage_2d_polar_pade)
   
   SET(ARGS ${CMAKE_BINARY_DIR}/gcsim2d_cartesian_input)
   ADD_TEST(NAME sim2d_gc_cart COMMAND test_2d_gc_cartesian ${ARGS})
   SET_TESTS_PROPERTIES(sim2d_gc_cart PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   
   SET(ARGS ${CMAKE_BINARY_DIR}/gcsim2d_polar_input)
   ADD_TEST(NAME sim2d_gc_polar COMMAND test_2d_gc_polar ${ARGS})
   SET_TESTS_PROPERTIES(sim2d_gc_polar PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET(ARGS ${CMAKE_BINARY_DIR}/vpsim2d_no_split_beam)
   ADD_TEST(NAME sim2d_vp_no_split COMMAND test_2d_vp_no_split ${ARGS})
   SET_TESTS_PROPERTIES(sim2d_vp_no_split PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET_TESTS_PROPERTIES(coordinate_transformations PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(fields_2d_alternative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(fields_1d_alternative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(characteristics_2d_explicit_euler PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(characteristics_1d_explicit_euler_conservative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(characteristics_2d_explicit_euler PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(characteristics_1d_trapezoid PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(characteristics_1d_trapezoid_conservative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(characteristics_2d_explicit_euler_conservative PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(characteristics_2d_verlet PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(advection_1d_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(
     advection_1d_non_uniform_cubic_splines
     PROPERTIES 
     PASS_REGULAR_EXPRESSION
     "PASSED")
   SET_TESTS_PROPERTIES(advection_2d_BSL PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(gyroaverage_polar_hermite PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(gyroaverage_polar_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   SET_TESTS_PROPERTIES(gyroaverage_polar_pade PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")



ADD_TEST(NAME electric_field_accumulators COMMAND test_e_field_accumulator_2d)

ADD_TEST(NAME mapped_meshes COMMAND test_mapped_meshes_1d
				    test_mapped_meshes_2d)

ADD_TEST(NAME ode_solvers COMMAND test_implicit_ode_nonuniform)

ADD_TEST(NAME BSL COMMAND bsl_1d_cubic_uniform_periodic
                          bsl_1d_cubic_nonuniform_periodic
                          bsl_1d_cubic_uniform_compact
                          bsl_1d_cubic_nonuniform_compact)
SET_TESTS_PROPERTIES(BSL PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


IF(FORTRANCL_FOUND)
   ADD_TEST(NAME opencl COMMAND test_opencl)
   SET_TESTS_PROPERTIES(opencl PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF(FORTRANCL_FOUND)
