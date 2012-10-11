PROJECT(slv2d)
ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )

SET(MPI_ENABLED           ON  CACHE BOOL  "Use parallel version of slv2d   ")
SET(FFTPACK_ENABLED       ON  CACHE BOOL  "Use dfftpack library            ")
SET(SILO_ENABLED          OFF CACHE BOOL  "Use SILO Library from LLNL      ")

FIND_PACKAGE(LAPACK  REQUIRED)

IF(FFTPACK_ENABLED)
   FIND_PACKAGE(FFTPACK)
ENDIF()

IF(FFTPACK_FOUND)
   ADD_DEFINITIONS(-D_FFTPACK)
ELSE()
   ADD_LIBRARY(sine vsint.f vsinti.f)
   ADD_DEFINITIONS(-D_SINE)
   SET(FFTPACK_LIBRARIES sine)
ENDIF()

IF(SILO_ENABLED) 
   FIND_PACKAGE(SILO  REQUIRED)
   IF(SILO_FOUND)
      ADD_DEFINITIONS(-D_SILO)
      INCLUDE_DIRECTORIES(${SILO_INCLUDE_DIRS})
      SET(CMAKE_EXE_LINKER_FLAGS "-lsupc++")
   ENDIF()
ENDIF()

IF(MPI_ENABLED)
   FIND_PACKAGE(MPI REQUIRED Fortran)
ENDIF()

IF(MPI_FOUND)
   MESSAGE(STATUS "MPI FOUND")
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   ADD_DEFINITIONS(-D_MPI)
ELSE(MPI_FOUND)
   MESSAGE(STATUS "MPI NOT FOUND, BUILD SEQUENTIAL VERSION")
ENDIF(MPI_FOUND)

MESSAGE( STATUS "PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR} )
MESSAGE( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )

SET(SELALIB_DIR ${PROJECT_SOURCE_DIR}/../../prototype/src)
MESSAGE( STATUS "SELALIB_DIR: " ${SELALIB_DIR} )

SET(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics")

FIND_PACKAGE(HDF5 REQUIRED Fortran)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})

SET( SLL_COMPONENTS precision 
                    utilities 
                    assert 
                    memory 
                    constants 
                    low_level_file_io)

FOREACH( component ${SLL_COMPONENTS} )
   MESSAGE(STATUS "${component}")
   INCLUDE_DIRECTORIES(${SELALIB_DIR}/${component})
   ADD_SUBDIRECTORY(${SELALIB_DIR}/${component} ${PROJECT_BINARY_DIR}/${component})
ENDFOREACH()


SET(CMAKE_Fortran_FLAGS "-ffree-line-length-none -fall-intrinsics")

SET(SLV2D_SRC silo_module.F90
              diagnostiquesm.F90
              diagnostiquesp.F90
              fft_module.F90
              geometry.F90
              maxwell2dfdtd.F90
              module_mpi.F90
              poisson2ddir.F90
              poisson2dpp.F90
              precision.F90
              splinenn.F90
              splinepp.F90
              splinepx.F90
              splinepy.F90
              transpose.F90
              vlasov1d.F90
              vlasov2d.F90 )

ADD_LIBRARY(slv2d STATIC ${SLV2D_SRC})

SET(LIBS sll_low_level_file_io
         sll_utilities
         sll_constants
         sll_memory
         sll_assertion
         sll_working_precision
         ${SILO_LIBRARIES} ${FFTPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES})

ADD_EXECUTABLE(vp2d vp2dinit.F90 vp2d.F90)
TARGET_LINK_LIBRARIES(vp2d slv2d ${LIBS} )
ADD_EXECUTABLE(vm2d vm2dinit.F90 vm2d.F90)
TARGET_LINK_LIBRARIES(vm2d slv2d ${LIBS} )

ENABLE_TESTING()

SET(ARGS ${CMAKE_SOURCE_DIR}/../input/slv2d.dat)
IF(MPI_FOUND)
   SET(PROCS 4)
   ADD_TEST(NAME test_vp2d COMMAND 
            ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
            ./vp2d
            ${MPIEXEC_POSTFLAGS} ${ARGS})
   ADD_TEST(NAME test_vm2d COMMAND 
            ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
            ./vm2d
            ${MPIEXEC_POSTFLAGS} ${ARGS})
ELSE()
   ADD_TEST(NAME vm2d COMMAND ./vm2d ${ARGS})
   ADD_TEST(NAME vp2d COMMAND ./vp2d ${ARGS})
ENDIF()

SET_TESTS_PROPERTIES(test_vp2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(test_vm2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

