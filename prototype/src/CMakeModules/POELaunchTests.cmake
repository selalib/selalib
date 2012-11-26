
IF(CMAKE_BUILD_TYPE STREQUAL DEBUG)
   ADD_TEST(NAME assert COMMAND test_assert)
   SET(passRegex "Assertion error triggered in file")
   SET_TESTS_PROPERTIES(assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}")
ENDIF(CMAKE_BUILD_TYPE STREQUAL DEBUG)

ADD_TEST(NAME memory COMMAND test_memory)
SET_TESTS_PROPERTIES(memory PROPERTIES TIMEOUT 20)
ADD_TEST(NAME constants COMMAND test_constants)
ADD_TEST(NAME utilities COMMAND test_utilities)
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
ADD_TEST(NAME integration COMMAND test_integration)

IF(FFTPACK_ENABLED)
  ADD_TEST(NAME periodic_interp COMMAND test_periodic_interp)
ENDIF()

ADD_TEST(NAME fft COMMAND test_fft)

IF(NOT STDF95)
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

ENDIF()

ADD_TEST(NAME electric_field_accumulators COMMAND test_e_field_accumulator_2d)

ADD_TEST(NAME WENO COMMAND test_WENO_interp test_WENO_recon)

ADD_TEST(NAME interpolators COMMAND test_interpolators_1d
				    test_interpolators_2d)

ADD_TEST(NAME quintic_interpolators_1d COMMAND test_quintic_interpolators_1d)

SET_TESTS_PROPERTIES(quintic_interpolators_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME odd_degree_interpolators_1d COMMAND test_odd_degree_interpolators_1d)

SET_TESTS_PROPERTIES(odd_degree_interpolators_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME mapped_meshes COMMAND test_mapped_meshes_1d
				    test_mapped_meshes_2d)

ADD_TEST(NAME fields COMMAND test_scalar_field)

ADD_TEST(NAME time_splitting COMMAND test_time_splitting)

ADD_TEST(NAME ode_solvers COMMAND test_implicit_ode_nonuniform)

ADD_TEST(NAME distribution_function COMMAND test_distribution_function)

ADD_TEST(NAME advection_field COMMAND test_advection_field)

ADD_TEST(NAME BSL COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_bsl_1d)
SET_TESTS_PROPERTIES(BSL PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME low_level_file_io COMMAND test_io)
SET_TESTS_PROPERTIES(low_level_file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME maxwell_2d_fdtd 
         COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_maxwell_2d_fdtd)
SET_TESTS_PROPERTIES(maxwell_2d_fdtd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

IF(FFTW_FOUND)
   ADD_TEST(NAME maxwell_2d_pstd 
            COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_maxwell_2d_pstd)
   SET_TESTS_PROPERTIES(maxwell_2d_pstd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()

#PARALLEL TESTS
IF(MPI_MODULE_ENABLED)

   IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)
      SET(PROCS 1)
      ADD_TEST(NAME poisson_periodic_cartesian_par_2d
               COMMAND
	       ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	       ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_poisson_2d_per_cart_par
	       ${MPIEXEC_POSTFLAGS} ${ARGS})
   ENDIF()

   SET(PROCS 16)
   ADD_TEST(NAME poisson_3d_periodic_par
         COMMAND
	 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_poisson_3d_periodic_par
	 ${MPIEXEC_POSTFLAGS} ${ARGS})

   SET_TESTS_PROPERTIES(poisson_3d_periodic_par 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   IF(NOT STDF95)

      SET(PROCS 16)
      ADD_TEST(NAME qns2d_with_finite_diff_par
            COMMAND
	    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_qns2d_with_finite_diff_par
	    ${MPIEXEC_POSTFLAGS} ${ARGS})
      SET_TESTS_PROPERTIES(qns2d_with_finite_diff_par 
                           PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

      SET(PROCS 16)
      ADD_TEST(NAME qns2d_angular_spectral_method_par
            COMMAND
	    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_qns2d_angular_spectral_method_par
	    ${MPIEXEC_POSTFLAGS} ${ARGS})
      SET_TESTS_PROPERTIES(qns2d_angular_spectral_method_par 
                           PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ENDIF()

   SET(PROCS 2)
   SET(ARGS "")
   ADD_TEST(NAME collective
            COMMAND
	   ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_collective
	   ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(collective PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS")

   SET(PROCS 8)
   SET(ARGS "")
   ADD_TEST(NAME remap_2d
	   COMMAND
	   ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_2d
	   ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET(PROCS 8)
   SET(ARGS "")
   ADD_TEST(NAME remap_3d
	   COMMAND
	   ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_3d
	   ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET(PROCS 16)
   ADD_TEST(NAME remap_4d
	   COMMAND
	   ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_4d
	   ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET(PROCS 8)
   ADD_TEST(NAME remap_6d
	   COMMAND
	   ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	   ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_6d
	   ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_6d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   IF(HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL)

      SET(PROCS 4)
      SET(ARGS "")
      ADD_TEST(NAME low_level_file_io_parallel COMMAND 
               ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_io_parallel
               ${MPIEXEC_POSTFLAGS} ${ARGS})
      SET_TESTS_PROPERTIES(low_level_file_io_parallel 
                           PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

      SET(PROCS 8)
      ADD_TEST(NAME test_vp4d_sim
            COMMAND
	    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vp4d_sim
	    ${MPIEXEC_POSTFLAGS} ${ARGS})
      SET_TESTS_PROPERTIES(test_vp4d_sim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ENDIF()

   IF(PASTIX_FOUND AND PTSCOTCH_FOUND AND MURGE_FOUND AND SCOTCH_FOUND)
      SET(PROCS 4)
      SET(ARGS 1000 3)
      ADD_TEST(NAME pastix 
               COMMAND
               ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_pastix ${ARGS})
      SET(pastix PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ENDIF()

ENDIF()

IF(FORTRANCL_FOUND)
  ADD_TEST(NAME opencl COMMAND test_opencl)
  SET_TESTS_PROPERTIES(opencl PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF(FORTRANCL_FOUND)

IF(FISHPACK_ENABLED)
   ADD_TEST(NAME fishpack4 COMMAND tpois3d)
   ADD_TEST(NAME fishpack3 COMMAND thwsssp)
   ADD_TEST(NAME fishpack2 COMMAND thwsplr)
   ADD_TEST(NAME fishpack1 COMMAND thwscyl)
ENDIF()

IF(MUDPACK_ENABLED AND Fortran_COMPILER STREQUAL "GFORTRAN")
   IF(MPI_MODULE_ENABLED AND HDF5_PARALLEL_ENABLED)
      SET(PROCS 4)
      SET(ARGS "")
      ADD_TEST(NAME multigrid COMMAND 
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
      ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_multigrid
      ${MPIEXEC_POSTFLAGS} ${ARGS})
      SET(multigrid PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ELSE()
      ADD_TEST(NAME mudpack COMMAND tmud34sp tmud24sp )
      SET(mudpack PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ENDIF()
ENDIF()
