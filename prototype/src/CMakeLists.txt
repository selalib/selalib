project(sll)
enable_language(Fortran)
cmake_minimum_required(VERSION 2.8)

# Specify the output for binary
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add new directory to find the FIND<PACKAGE>.cmake files. CMAKE_SOURCE_DIR is
# the directory where the project sources are located. CMAKE_MODULE_PATH is 
# the list of directories in which to search for CMake modules.
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )

# make sure that the default is a RELEASE. Here is where we need to set
# the parameters for the DEBUG or RELEASE versions. For instance, DEBUG should
# turn on assertions.
if (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Define folder to be deleted by make clean
#set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/modules/)
#set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/bin/)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/include/)
add_custom_target( clear
                   COMMAND find . -name *.mod -delete
		   COMMAND find . -name *.a -delete
		   COMMAND find . -name *.o -delete
		   COMMAND find . -name *.cmake -delete
		   COMMAND rm -rf bin
		   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		   COMMENT "Delete all .a, .mod, .o and folders /modules, /bin" VERBATIM)

INCLUDE(${CMAKE_SOURCE_DIR}/CMakeModules/FortranConfig.cmake)

INCLUDE(${CMAKE_SOURCE_DIR}/CMakeModules/HDF5Config.cmake)

INCLUDE(${CMAKE_SOURCE_DIR}/CMakeModules/FFTConfig.cmake)

INCLUDE(${CMAKE_SOURCE_DIR}/CMakeModules/MPIConfig.cmake)

INCLUDE(${CMAKE_SOURCE_DIR}/CMakeModules/DoxygenConfig.cmake)

INCLUDE(${CMAKE_SOURCE_DIR}/CMakeModules/LapackConfig.cmake)


# include's .h files
INCLUDE_DIRECTORIES(remap 
                    poisson_solvers 
                    utilities 
                    memory 
                    assert 
                    splines 
                    numerical_utilities 
                    precision 
                    fields 
                    fft 
                    remap) 



ENABLE_TESTING()

######################################################################
#  III. SUBDIRECTORY                                                 #
######################################################################

#low-level
ADD_SUBDIRECTORY(memory)
ADD_TEST(NAME memory COMMAND test_memory)
SET_TESTS_PROPERTIES(memory PROPERTIES TIMEOUT 20)

##########################################################

ADD_SUBDIRECTORY(precision)
#ADD_TEST(precision precision/test_precision) #too long

##########################################################

ADD_SUBDIRECTORY(assert)
ADD_TEST(NAME assert COMMAND test_assert)
SET(passRegex "Assertion error triggered in file")
SET_TESTS_PROPERTIES(assert PROPERTIES PASS_REGULAR_EXPRESSION "${passRegex}")

##########################################################

ADD_SUBDIRECTORY(constants)
ADD_TEST(NAME constants COMMAND test_constants)

##########################################################

ADD_SUBDIRECTORY(utilities)
ADD_TEST(NAME utilities COMMAND test_utilities)

##########################################################

ADD_SUBDIRECTORY(low_level_file_io)
ADD_TEST(NAME low_level_file_io COMMAND test_io)
SET_TESTS_PROPERTIES(low_level_file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
IF(HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)
   SET(PROCS 4)
   SET(ARGS "")
   ADD_TEST(NAME low_level_file_io_parallel COMMAND 
            ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_io_parallel
            ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(low_level_file_io_parallel 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()

##########################################################

ADD_SUBDIRECTORY(timer)
ADD_TEST(NAME timer COMMAND test_timer)

##########################################################

#mid-level
ADD_SUBDIRECTORY(numerical_utilities)
ADD_TEST(NAME tridiagonal COMMAND test_tridiagonal)
ADD_TEST(NAME lagrange COMMAND test_lagrange)

ADD_SUBDIRECTORY(non_linear_solvers)
ADD_TEST(NAME toeplitz_penta_diagonal COMMAND test_toeplitz_penta_diagonal)
SET_TESTS_PROPERTIES(toeplitz_penta_diagonal PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME newton_raphson COMMAND test_newton_raphson)

##########################################################

ADD_SUBDIRECTORY(splines)
ADD_TEST(NAME splines COMMAND test_splines ) 
SET_TESTS_PROPERTIES(splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME splines_arbitrary_degree COMMAND test_arbitrary_degree_splines)
SET_TESTS_PROPERTIES(splines_arbitrary_degree PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
# For Michel: add here the test_non_unif_splines... 
ADD_TEST(NAME quintic_splines COMMAND test_quintic_splines)
SET_TESTS_PROPERTIES(quintic_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ADD_TEST(NAME odd_degree_splines COMMAND test_odd_degree_splines)
SET_TESTS_PROPERTIES(odd_degree_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

##########################################################

ADD_SUBDIRECTORY(WENO)
ADD_SUBDIRECTORY(integration)
ADD_TEST(NAME integration COMMAND test_integration)

##########################################################

IF(FFTPACK_ENABLED)
   ADD_SUBDIRECTORY(fftpack)
ENDIF()

##########################################################

ADD_SUBDIRECTORY(fft)
ADD_TEST(NAME fft COMMAND test_fft)

##########################################################

IF(MPI_MODULE_ENABLED)
   ADD_SUBDIRECTORY(collective)
   ADD_SUBDIRECTORY(remap)
   # we can call ${MPIEXEC_MAX_NUMPROCS} to have the number of procs
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
   ADD_TEST(NAME remap
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap
	${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ADD_TEST(NAME remap_2d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_2d
	${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   # to make the test meaningful at higher dimensions, we need more processes...
   SET(PROCS 16)
   ADD_TEST(NAME remap_4d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_4d
	${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET(PROCS 16)
   ADD_TEST(NAME remap_6d
	COMMAND
	${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_remap_6d
	${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_6d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET(PROCS 16)
   ADD_TEST(NAME poisson_periodic_cartesian_par_2d
         COMMAND
	 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_poisson_2d_per_cart_par
	 ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_6d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()

##########################################################

IF(NOT STDF95)
   ADD_SUBDIRECTORY(poisson_solvers)
   IF(FFTPACK_ENABLED)
      ADD_TEST(NAME poisson_solvers COMMAND test_poisson_1d)
   ENDIF(FFTPACK_ENABLED)

#ADD_TEST(poisson_solvers poisson_solvers/test_poisson_2d)
#ADD_TEST(poisson_solvers poisson_solvers/test_poisson_3d_seq)
#ADD_TEST(poisson_solvers mpirun -np 2 poisson_solvers/test_poisson_3d_par)
#ADD_TEST(poisson_solvers poisson_solvers/unit_test)

   ADD_SUBDIRECTORY(quasi_neutral_solvers)
   ADD_TEST(NAME qns2d_with_finite_diff_seq 
            COMMAND test_qns2d_with_finite_diff_seq)
   SET_TESTS_PROPERTIES(qns2d_with_finite_diff_seq 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ADD_TEST(NAME qns2d_angular_spectral_method_seq 
            COMMAND test_qns2d_angular_spectral_method_seq)
   SET_TESTS_PROPERTIES(qns2d_angular_spectral_method_seq 
                        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
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

ADD_SUBDIRECTORY(electric_field_accumulators)
ADD_TEST(NAME electric_field_accumulators COMMAND test_e_field_accumulator_2d)

##########################################################
#high-level
ADD_SUBDIRECTORY(interpolators)
ADD_TEST(NAME WENO COMMAND test_WENO_interp test_WENO_recon)
ADD_TEST(NAME interpolators COMMAND test_interpolators_1d
				    test_interpolators_2d)

ADD_TEST(NAME quintic_interpolators_1d COMMAND test_quintic_interpolators_1d)
SET_TESTS_PROPERTIES(quintic_interpolators_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_TEST(NAME odd_degree_interpolators_1d COMMAND test_odd_degree_interpolators_1d)
SET_TESTS_PROPERTIES(odd_degree_interpolators_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

##########################################################

ADD_SUBDIRECTORY(mapped_meshes)
ADD_TEST(NAME mapped_meshes COMMAND test_mapped_meshes_1d
				    test_mapped_meshes_2d)

##########################################################

ADD_SUBDIRECTORY(fields)
ADD_TEST(NAME fields COMMAND test_scalar_field)

##########################################################

ADD_SUBDIRECTORY(time_solvers)
ADD_TEST(NAME time_splitting COMMAND test_time_splitting)

##########################################################

ADD_SUBDIRECTORY(ode_solvers)
ADD_TEST(NAME ode_solvers COMMAND test_implicit_ode_nonuniform)

##########################################################

ADD_SUBDIRECTORY(distribution_function)
ADD_TEST(NAME distribution_function COMMAND test_distribution_function)

##########################################################

ADD_SUBDIRECTORY(advection_field)
ADD_TEST(NAME advection_field COMMAND test_advection_field)

##########################################################

ADD_SUBDIRECTORY(LinRood)

##########################################################

ADD_SUBDIRECTORY(maxwell_solvers)
ADD_TEST(NAME maxwell_2d_fdtd 
         COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_maxwell_2d_fdtd)
IF(FFTW_FOUND)
   ADD_TEST(NAME maxwell_2d_pstd 
            COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_maxwell_2d_pstd)
ENDIF()
SET_TESTS_PROPERTIES(maxwell_2d_fdtd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(maxwell_2d_pstd PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

##########################################################

ADD_SUBDIRECTORY(BSL)
ADD_TEST(NAME BSL COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_bsl_1d)
SET_TESTS_PROPERTIES(BSL PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
#ADD_TEST(NAME landau_2d COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_bsl_landau_2d)
#SET_TESTS_PROPERTIES(landau_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


##########################################################

ADD_SUBDIRECTORY(CSL)
#ADD_TEST(CSL CSL/test_csl)
#ADD_TEST(simulation_1d simulation/test_simulation_1d)

#ADD_SUBDIRECTORY(quasi_neutral_solvers) don't compile problem with fftw3_thread! top
#ADD_TEST(quasi_neutral_solvers quasi_neutral_solvers/test_quasi_neutral_solver)

##########################################################

ADD_SUBDIRECTORY(selalib_level EXCLUDE_FROM_ALL)
#add_custom_target(install make selalib_low_level)

##########################################################

ADD_SUBDIRECTORY(simulation)
IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)
  SET(PROCS 8)
  ADD_TEST(NAME test_vp4d_sim
         COMMAND
	 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vp4d_sim
	 ${MPIEXEC_POSTFLAGS} ${ARGS})
  SET_TESTS_PROPERTIES(test_vp4d_sim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()

##########################################################

SET(FORTRANCL_ENABLED NO CACHE BOOL "Enable OpenCL test programs")
# recherche de FORTRANCL
IF(FORTRANCL_ENABLED)
  FIND_PACKAGE(FORTRANCL QUIET)
  IF(FORTRANCL_FOUND)
    MESSAGE(STATUS "FORTRANCL FOUND")
    INCLUDE_DIRECTORIES(${FORTRANCL_INCLUDE_DIRS})
    ADD_SUBDIRECTORY(opencl)
    ADD_TEST(NAME opencl COMMAND test_opencl)
    SET_TESTS_PROPERTIES(opencl PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  ENDIF(FORTRANCL_FOUND)
ENDIF(FORTRANCL_ENABLED)

##########################################################

SET(PASTIX_ENABLED NO CACHE BOOL " ")

IF(PASTIX_ENABLED)
   IF( DEFINED ENV{SCOTCH_ROOT} )
      SET(SCOTCH_ROOT $ENV{SCOTCH_ROOT})
      SET(PTSCOTCH_ROOT $ENV{SCOTCH_ROOT})
   ENDIF()
   IF( DEFINED ENV{PASTIX_ROOT} )
      SET(PASTIX_ROOT $ENV{PASTIX_ROOT})
      SET(MURGE_ROOT $ENV{PASTIX_ROOT})
   ENDIF()
   FIND_PACKAGE(RT QUIET)
   FIND_PACKAGE(SCOTCH QUIET)
   FIND_PACKAGE(PTSCOTCH QUIET)
   FIND_PACKAGE(MURGE QUIET)
   FIND_PACKAGE(PASTIX QUIET)
   IF(PASTIX_FOUND AND PTSCOTCH_FOUND AND MURGE_FOUND AND SCOTCH_FOUND)
      message(STATUS "PASTIX FOUND")
      INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})
      ADD_SUBDIRECTORY(pastix)
      SET(PROCS 4)
      SET(ARGS 1000 3)
      ADD_TEST(NAME pastix 
               COMMAND
               ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_pastix ${ARGS})
      SET(pastix PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ELSE()
      message(STATUS "PASTIX NOT FOUND")
   ENDIF()
ENDIF(PASTIX_ENABLED)

##########################################################

SET(FISHPACK_ENABLED OFF CACHE BOOL " ")
IF(FISHPACK_ENABLED)
   ADD_SUBDIRECTORY(fishpack)
   ADD_TEST(NAME fishpack4 COMMAND tpois3d)
   ADD_TEST(NAME fishpack3 COMMAND thwsssp)
   ADD_TEST(NAME fishpack2 COMMAND thwsplr)
   ADD_TEST(NAME fishpack1 COMMAND thwscyl)
ENDIF()

##########################################################

SET(MUDPACK_ENABLED OFF CACHE BOOL "Enable multigrid solvers")
IF(MUDPACK_ENABLED AND Fortran_COMPILER STREQUAL "GFORTRAN")
   IF(MPI_MODULE_ENABLED AND HDF5_PARALLEL_ENABLED)
      ADD_SUBDIRECTORY(multigrid)
   ELSE()
      ADD_SUBDIRECTORY(mudpack)
   ENDIF()
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
