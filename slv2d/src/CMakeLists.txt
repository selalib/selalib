project(slv2d)
enable_language(Fortran)
cmake_minimum_required(VERSION 2.8)

set(HDF5_ENABLED ON CACHE BOOL "Use HDF5 format for data output ")
set(HDF5_PARALLEL_ENABLED OFF CACHE BOOL "Use Parallel HDF5")

IF(HDF5_ENABLED)
	IF( DEFINED ENV{HDF5_ROOT} )
		SET(HDF5_ROOT $ENV{HDF5_ROOT})
	ENDIF()
	FIND_PACKAGE(HDF5 QUIET)
	IF(HDF5_FOUND)
		MESSAGE(STATUS "HDF5 FOUND")
		INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
	ENDIF()
ENDIF(HDF5_ENABLED)

# add the cache entry FFTPACK_ENABLED for enable/disable fftpack
SET(FFTPACK_ENABLED OFF CACHE BOOL " ")

FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(ZLIB REQUIRED)

#add flag to disable MPI stuff for debug
SET(MPI_MODULE_ENABLED ON CACHE BOOL " ")

FIND_PACKAGE(MPI)
IF(MPI_FOUND)
   MESSAGE(STATUS "MPI FOUND")
   INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   SET(CMAKE_Fortran_FLAGS "-D_MPI")
ELSE(MPI_FOUND)
   MESSAGE(STATUS "MPI NOT FOUND")
   SET(MPI_MODULE_ENABLED OFF CACHE BOOL " " FORCE)
ENDIF(MPI_FOUND)

FIND_LIBRARY(FFTPACK_LIBRARIES
		NAMES dfftpack
		HINTS /usr/local
		PATH_SUFFIXES lib 
		DOC "PATH TO libdfftpack.a")

FIND_PATH(SILO_INCLUDE_DIRS
		NAMES silo.inc
		HINTS /usr/local
		PATH_SUFFIXES include 
		DOC "PATH TO silo.inc")

FIND_LIBRARY(SILO_LIBRARIES
		NAMES silo siloh5
		HINTS /usr/local
		PATH_SUFFIXES lib 
		DOC "PATH TO libsilo.a or libsiloh5.a")

if (SILO_INCLUDE_DIRS AND SILO_LIBRARIES)
  set(SILO_FOUND YES)
endif(SILO_INCLUDE_DIRS AND SILO_LIBRARIES)

SET(SLV2D_SRC diagnostiquesm.F90
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

ADD_EXECUTABLE(vp2d vp2dinit.F90 vp2d.F90)
TARGET_LINK_LIBRARIES(vp2d slv2d ${SILO_LIBRARIES} ${FFTPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} -lsupc++)
ADD_EXECUTABLE(vm2d vm2dinit.F90 vm2d.F90)
TARGET_LINK_LIBRARIES(vm2d slv2d ${SILO_LIBRARIES} ${FFTPACK_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} -lsupc++)

ENABLE_TESTING()

ADD_TEST(NAME vm2d COMMAND test_vm2d)

IF(HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)
   SET(PROCS 4)
   SET(ARGS "")
   ADD_TEST(NAME vp2d COMMAND 
            ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vp2d
            ${MPIEXEC_POSTFLAGS} ${ARGS})
   SET_TESTS_PROPERTIES(vp2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()
