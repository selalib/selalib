PROJECT(slv2d)
ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )
SET(SELALIB_DIR ${PROJECT_SOURCE_DIR}/../../prototype/src)
MESSAGE( STATUS "SELALIB_DIR: " ${SELALIB_DIR} )
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules 
                       ${SELALIB_DIR}/CMakeModules
                       ${CMAKE_MODULE_PATH} )

INCLUDE(${CMAKE_SOURCE_DIR}/CMakeModules/FortranConfig.cmake)

SET(FFTPACK_ENABLED ON  CACHE BOOL "Use dfftpack library       ")
SET(SILO_ENABLED    OFF CACHE BOOL "Use SILO Library from LLNL ")
SET(FFTW_ENABLED    OFF CACHE BOOL "Use fftw library       ")

INCLUDE(${SELALIB_DIR}/CMakeModules/LapackConfig.cmake)

IF(SILO_ENABLED) 
   FIND_PACKAGE(SILO  REQUIRED)
   IF(SILO_FOUND)
      ADD_DEFINITIONS(-D_SILO)
      INCLUDE_DIRECTORIES(${SILO_INCLUDE_DIRS})
      SET(CMAKE_EXE_LINKER_FLAGS "-lsupc++")
   ENDIF()
ENDIF()

FIND_PACKAGE(MPI REQUIRED Fortran)
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

SET(SELALIB_DIR ${PROJECT_SOURCE_DIR}/../../prototype/src)
MESSAGE( STATUS "SELALIB_DIR: " ${SELALIB_DIR} )

IF(FFTPACK_ENABLED)
   FIND_PACKAGE(FFTPACK)
   IF(FFTPACK_FOUND)
      ADD_DEFINITIONS(-D_FFTPACK)
   ENDIF()
ENDIF()

IF(FFTW_ENABLED)
   FIND_PACKAGE(FFTW)
   IF(FFTW_FOUND)
      ADD_DEFINITIONS(-D_FFTW)
      INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
   ENDIF()
ENDIF()

#------------------------------------
# Check if using IBM compilers
#------------------------------------
IF($ENV{FC} MATCHES "xlf")
  INCLUDE(IBMCompilers)
ENDIF()

SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "Enable HDF5 parallel")

SET(HDF5_ROOT $ENV{HDF5_ROOT})
FIND_PACKAGE(HDF5 REQUIRED Fortran)
SET (HDF5_ENABLE_PARALLEL @HDF5_ENABLE_PARALLEL@)
IF(HDF5_ENABLE_PARALLEL) 
   MESSAGE(STATUS "HDF5 parallel supported")
ELSE()
   MESSAGE(SEND_ERROR "HDF5 parallel not supported")
ENDIF()

IF(HDF5_INCLUDE_DIRS)
   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS}/fortran)
ENDIF()
IF(HDF5_INCLUDE_DIR_FORTRAN)
   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR_FORTRAN})
ENDIF()

#Pas trouver de meilleur moyen pour l'instant...
FILE(COPY ${SELALIB_DIR}/fields/geometry_functions.F90 DESTINATION ${PROJECT_SOURCE_DIR}/fields/)

SET( SLL_COMPONENTS precision 
                    utilities 
                    assert 
                    memory 
                    constants 
                    timer
                    numerical_utilities
                    electric_field_accumulators
                    remap
                    splines
                    WENO
                    interpolators
                    collective
                    low_level_file_io)

FOREACH( component ${SLL_COMPONENTS} )
   MESSAGE(STATUS "${component}")
   INCLUDE_DIRECTORIES(${SELALIB_DIR}/${component})
   ADD_SUBDIRECTORY(${SELALIB_DIR}/${component} ${PROJECT_BINARY_DIR}/${component})
ENDFOREACH()

IF(Fortran_COMPILER STREQUAL "GFORTRAN")
   SET(CMAKE_Fortran_FLAGS "-g -Wall -ffree-line-length-none -fall-intrinsics -fbounds-check")
ENDIF()

IF(Fortran_COMPILER STREQUAL "XLF")
   SET(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ENDIF()

IF(Fortran_COMPILER STREQUAL "IFORT")

ENDIF()

SET(SLV2D_SRC silo_module.F90
              diagnostiques.F90
              fft_module.F90
              geometry.F90
              maxwell2dfdtd.F90
              poisson2dpp_seq.F90
              poisson2dpp_par.F90
              poisson2d_periodic.F90
              precision.F90
              splinenn.F90
              splinepp.F90
              splinepx.F90
              splinepy.F90
              transpose.F90
              vlasov1d.F90
              vlasov2d.F90 )

ADD_LIBRARY(slv2d STATIC ${SLV2D_SRC})

ADD_DEPENDENCIES(slv2d sll_precision
                       sll_assertion
                       sll_constants
                       sll_memory
                       sll_low_level_file_io
                       sll_timer
                       sll_collective
                       sll_remap
                       sll_splines
                       sll_interpolators
                       sll_numerical_utilities
                       sll_low_level_file_io
                       sll_low_level_file_io_parallel)

SET(LIBS sll_utilities
         sll_constants
         sll_memory
         sll_assertion
         sll_working_precision
         sll_collective
         sll_interpolators
         sll_splines
         sll_remap
         sll_numerical_utilities
         sll_low_level_file_io
         sll_low_level_file_io_parallel
         ${SILO_LIBRARIES} 
         ${LAPACK_LIBRARIES} 
         ${FFTPACK_LIBRARIES} 
         ${FFTW_LIBRARIES}
         ${HDF5_LIBRARIES})

ADD_EXECUTABLE(vp2d vp2dinit.F90 vp2d.F90)
TARGET_LINK_LIBRARIES(vp2d slv2d ${LIBS})

ADD_EXECUTABLE(vm2d vm2dinit.F90 vm2d.F90)
TARGET_LINK_LIBRARIES(vm2d slv2d ${LIBS})

ADD_EXECUTABLE(test_poisson2d_par test_poisson2d_par.F90)
TARGET_LINK_LIBRARIES(test_poisson2d_par slv2d ${LIBS})

ADD_EXECUTABLE(vp2d_cubic vp2d_selalib.F90 sll_vlasov2d.F90)
TARGET_LINK_LIBRARIES(vp2d_cubic slv2d ${LIBS})

ADD_EXECUTABLE(vp2d_keen vp2d_keen.F90 sll_vlasov2d.F90)
TARGET_LINK_LIBRARIES(vp2d_keen slv2d ${LIBS})

#FILE(REMOVE ${CMAKE_BINARY_DIR}/modules/sll_vlasov2d.mod)
#ADD_DEFINITIONS(-D_QUINTIC)
#ADD_EXECUTABLE(vp2d_quintic vp2d_selalib.F90 sll_vlasov2d.F90)
#TARGET_LINK_LIBRARIES(vp2d_quintic slv2d ${LIBS})

FILE(REMOVE ${CMAKE_BINARY_DIR}/modules/sll_vlasov2d.mod)
ADD_DEFINITIONS(-D_TWO_D)
ADD_EXECUTABLE(vp2d_sll_interp2d vp2d_selalib.F90 sll_vlasov2d.F90)
TARGET_LINK_LIBRARIES(vp2d_sll_interp2d slv2d ${LIBS})

IF (FFTW_FOUND) 
   ADD_EXECUTABLE(fftw_solver fftw_solver.F90)
   TARGET_LINK_LIBRARIES(fftw_solver ${LIBS})
   IF(FFTW_MPI_LIBRARY)
      ADD_EXECUTABLE(test_fftw_mpi test_fftw_mpi.F90 fftw3.F90)
      TARGET_LINK_LIBRARIES(test_fftw_mpi ${FFTW_LIBRARIES})
      ADD_EXECUTABLE(fftw_solver_mpi fftw_solver_mpi.F90 fftw3.F90)
      TARGET_LINK_LIBRARIES(fftw_solver_mpi ${FFTW_LIBRARIES} sll_memory)
   ENDIF()
ENDIF()

ENABLE_TESTING()

SET(PROCS 1)
ADD_TEST(NAME test_poisson2d_par COMMAND 
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
         ./test_poisson2d_par
         ${MPIEXEC_POSTFLAGS} ${ARGS})

SET(PROCS 4)
SET(ARGS ${CMAKE_SOURCE_DIR}/../input/slv2d.nml)

ADD_TEST(NAME test_vp2d COMMAND 
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
         ./vp2d
         ${MPIEXEC_POSTFLAGS} ${ARGS})

ADD_TEST(NAME test_vm2d COMMAND 
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
         ./vm2d
         ${MPIEXEC_POSTFLAGS} ${ARGS})

ADD_TEST(NAME test_vp2d_cubic COMMAND 
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
         ./vp2d_cubic
         ${MPIEXEC_POSTFLAGS} ${ARGS})

#ADD_TEST(NAME test_vp2d_quintic COMMAND 
#         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
#         ./vp2d_quintic
#         ${MPIEXEC_POSTFLAGS} ${ARGS})

SET_TESTS_PROPERTIES(test_vp2d          PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(test_vm2d          PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(test_poisson2d_par PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(test_vp2d_cubic    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
