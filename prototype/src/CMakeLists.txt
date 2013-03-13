IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} 
       CACHE STRING "Choose the type of build, options are: Debug Release.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug 
       CACHE STRING "Choose the type of build, options are: Debug Release.")
ENDIF()

IF(SLL_CONFIGURE_OPTIONS_FILE)
INCLUDE(${SLL_CONFIGURE_OPTIONS_FILE})
ENDIF()

PROJECT(SLL)

ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Specify the output for binary
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add new directory to find the FIND<PACKAGE>.cmake files. CMAKE_SOURCE_DIR is
# the directory where the project sources are located. CMAKE_MODULE_PATH is 
# the list of directories in which to search for CMake modules.
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )

IF(CMAKE_BUILD_TYPE MATCHES Debug)
   ADD_DEFINITIONS(-DDEBUG)
ENDIF()

# Define folder to be deleted by make clean
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/modules/)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/bin/)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/include/)

ADD_CUSTOM_TARGET( clear
   COMMAND find . -name *.mod -delete
   COMMAND find . -name *.a -delete
   COMMAND find . -name *.o -delete
   COMMAND find . -name *.cmake -delete
   COMMAND rm -rf bin
   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
   COMMENT "Delete all .a, .mod, .o and folders /modules, /bin" VERBATIM)

SET(STDF95_ENABLED OFF CACHE BOOL "Use F95 norm")

INCLUDE(FortranConfig)
INCLUDE(LapackConfig)
INCLUDE(HDF5Config)
INCLUDE(FFTConfig)
INCLUDE(MPIConfig)
INCLUDE(DoxygenConfig)

# include's .h files
INCLUDE_DIRECTORIES(utilities 
                    memory 
                    assert 
                    splines 
                    numerical_utilities 
                    precision 
                    fields 
                    fft) 
ENABLE_TESTING()

 ADD_SUBDIRECTORY(memory)
 ADD_SUBDIRECTORY(precision)
 ADD_SUBDIRECTORY(assert)
 ADD_SUBDIRECTORY(constants)
 ADD_SUBDIRECTORY(utilities)
 ADD_SUBDIRECTORY(low_level_file_io)
 ADD_SUBDIRECTORY(timer)
 ADD_SUBDIRECTORY(numerical_utilities)
 ADD_SUBDIRECTORY(non_linear_solvers)
 ADD_SUBDIRECTORY(splines)
 ADD_SUBDIRECTORY(WENO)
 ADD_SUBDIRECTORY(integration)
 ADD_SUBDIRECTORY(fft)
IF(FFTPACK_ENABLED)
   ADD_SUBDIRECTORY(fftpack)
   ADD_SUBDIRECTORY(periodic_interpolation)
ENDIF()
IF(MPI_MODULE_ENABLED)
   ADD_SUBDIRECTORY(collective)
   ADD_SUBDIRECTORY(remap)
ENDIF()
 ADD_SUBDIRECTORY(poisson_solvers)
 ADD_SUBDIRECTORY(quasi_neutral_solvers)
 ADD_SUBDIRECTORY(electric_field_accumulators)
 ADD_SUBDIRECTORY(interpolators)
 ADD_SUBDIRECTORY(mapped_meshes)
 ADD_SUBDIRECTORY(fields)
IF(NOT STDF95)
 ADD_SUBDIRECTORY(time_solvers)
ENDIF(NOT STDF95)
 ADD_SUBDIRECTORY(ode_solvers)
IF (NOT STDF95)
 ADD_SUBDIRECTORY(distribution_function)
 ADD_SUBDIRECTORY(advection_field)
 ADD_SUBDIRECTORY(LinRood)
ENDIF(NOT STDF95)
 ADD_SUBDIRECTORY(maxwell_solvers)
 ADD_SUBDIRECTORY(BSL)
IF(NOT STDF95)
 ADD_SUBDIRECTORY(CSL)
ENDIF(NOT STDF95)
 ADD_SUBDIRECTORY(selalib_level)
 ADD_SUBDIRECTORY(visu_pic)
 ADD_SUBDIRECTORY(simulation)

##########################################################

#ADD_SUBDIRECTORY(simulation)
#IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)
#  SET(PROCS 8)
#  ADD_TEST(NAME test_vp4d_sim
#         COMMAND
#	 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
#	 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vp4d_sim
#	 ${MPIEXEC_POSTFLAGS} ${ARGS})
#  SET_TESTS_PROPERTIES(test_vp4d_sim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
#  SET(PROCS 16)
#  ADD_TEST(NAME test_vp6d_sim
#         COMMAND
#	 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
#	 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vp6d_sim
#	 ${MPIEXEC_POSTFLAGS} ${ARGS})
#  SET_TESTS_PROPERTIES(test_vp6d_sim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

#ENDIF()

##########################################################

SET(FORTRANCL_ENABLED NO CACHE BOOL "Enable OpenCL test programs")

IF(FORTRANCL_ENABLED)
  FIND_PACKAGE(FORTRANCL QUIET)
  IF(FORTRANCL_FOUND)
    MESSAGE(STATUS "FORTRANCL FOUND")
    INCLUDE_DIRECTORIES(${FORTRANCL_INCLUDE_DIRS})
    ADD_SUBDIRECTORY(opencl)
  ENDIF(FORTRANCL_FOUND)
ENDIF(FORTRANCL_ENABLED)

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
      MESSAGE(STATUS "PASTIX FOUND")
      INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})
      ADD_SUBDIRECTORY(pastix)
   ELSE()
      MESSAGE(STATUS "PASTIX NOT FOUND")
   ENDIF()
ENDIF(PASTIX_ENABLED)

SET(FISHPACK_ENABLED OFF CACHE BOOL " ")

IF(FISHPACK_ENABLED)
   ADD_SUBDIRECTORY(fishpack)
ENDIF()

SET(MUDPACK_ENABLED OFF CACHE BOOL "Enable multigrid solvers")

IF(MUDPACK_ENABLED)
   IF(HDF5_PARALLEL_ENABLED)
      ADD_SUBDIRECTORY(multigrid)
   ELSE()
      ADD_SUBDIRECTORY(mudpack)
   ENDIF()
ENDIF()

MESSAGE(STATUS "############################################################")
MESSAGE(STATUS "####")
MESSAGE(STATUS "#### CMAKE_BUILD_TYPE:${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "#### CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
MESSAGE(STATUS "#### Fortran FLAGS ")
IF(CMAKE_BUILD_TYPE MATCHES "Debug") 
   MESSAGE(STATUS "#### ${CMAKE_Fortran_FLAGS_DEBUG}")
ELSE()
   MESSAGE(STATUS "#### ${CMAKE_Fortran_FLAGS_RELEASE}")
ENDIF()
MESSAGE(STATUS "####")
MESSAGE(STATUS "############################################################")

IF(STDF95)
MESSAGE(STATUS " ")
MESSAGE(STATUS "####################  F95  ###################")
MESSAGE(STATUS "## Directories not compiled in Fortran 95 : ##")
MESSAGE(STATUS "## advection_field                          ##")
MESSAGE(STATUS "## Linrood                                  ##")
MESSAGE(STATUS "## distribution_function                    ##")
MESSAGE(STATUS "## time_solvers                             ##")
MESSAGE(STATUS "## fields (except geometry_functions.F90)   ##")
MESSAGE(STATUS "##############################################")
ENDIF(STDF95)

INCLUDE(SequentialTests)


IF(MPI_MODULE_ENABLED)
   INCLUDE(ProcessorCount)
   INCLUDE(POEConfig)
   INCLUDE(ParallelTests)
ENDIF()

SET(PETSC_ENABLED OFF CACHE BOOL "Enable PETSc solvers")
IF(PETSC_ENABLED)
   FIND_PACKAGE(PETSc)
   IF(PETSC_FOUND)
      INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
      ADD_SUBDIRECTORY(petsc_solvers)
   ENDIF(PETSC_FOUND)
ENDIF()
