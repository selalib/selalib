IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: Debug Release.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: Debug Release.")
ENDIF()

PROJECT(GYSELA_PUBLIC)
ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

MESSAGE( STATUS "GYSELA_PUBLIC_BINARY_DIR: ${GYSELA_PUBLIC_BINARY_DIR}" )
MESSAGE( STATUS "GYSELA_PUBLIC_SOURCE_DIR: ${GYSELA_PUBLIC_SOURCE_DIR}" )
IF (CMAKE_PROJECT_NAME MATCHES "GYSELA_PUBLIC")
   SET(SLL_SOURCE_DIR "${GYSELA_PUBLIC_SOURCE_DIR}/../../prototype/src")
ENDIF()
MESSAGE( STATUS "SLL_SOURCE_DIR: ${SLL_SOURCE_DIR}" )

SET (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake 
                       ${SLL_SOURCE_DIR}/CMakeModules
                       ${CMAKE_MODULE_PATH} )

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

INCLUDE_DIRECTORIES(${SLL_SOURCE_DIR}/selalib_level)
LINK_DIRECTORIES(${SLL_SOURCE_DIR}/selalib_level)

MESSAGE(STATUS "CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
INCLUDE(F90FLAGS)
MESSAGE(STATUS "CMAKE_Fortran_FLAGS:${CMAKE_Fortran_FLAGS}")

IF(NOT MPI_FOUND)
   FIND_PACKAGE(MPI REQUIRED Fortran)
ENDIF(NOT MPI_FOUND)
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "Enable HDF5 parallel")

IF(NOT HDF5_FOUND)
   INCLUDE(HDF5Config)
   IF(HDF5_ENABLE_PARALLEL) 
      MESSAGE(STATUS "HDF5 parallel supported")
   ELSE()
      MESSAGE(SEND_ERROR "HDF5 parallel not supported")
      MESSAGE(SEND_ERROR "./configure --enable-fortran --enable-parallel CC=mpicc FC=mpif90")
   ENDIF()
ENDIF(NOT HDF5_FOUND)

IF (CMAKE_PROJECT_NAME MATCHES "GYSELA_PUBLIC")
   INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../../prototype/src/assert)
   INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../../prototype/src/memory)
   INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../../prototype/src/precision)
   INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/../../prototype/src/utilities)
   ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/../../prototype/src/selalib_level ${CMAKE_BINARY_DIR}/selalib_level EXCLUDE_FROM_ALL)
ELSE()
   INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
ENDIF()

SET(GYSELA_SRCS1 prec_const.f90 globals.f90 clock.f90 mem_alloc.f90 
	          utils.f90 MPIutils.f90 OMPutils.f90 
	          geometry.f90 gys_alloc.f90 local_spline.f90)

SET(GYSELA_SRCS2 HDF5_io.f90 resu3D_saving.f90 f5D_saving.f90)

SET(GYSELA_SRCS3 bessel.f90 LUsolving.f90 read_data.f90 spline1d.f90 
	          interpolation.f90 coord_system.f90 integration.f90
	          read_profile.f90 init_profile.f90 init_current.f90 
	          init_magnetic.f90 
	          fft_NRF90.f90 fft.f90 gyroaverage.f90 
	          fequil4d.f90 fdistribu5d.f90 
	          poisson.f90 efield.f90 
	          filter.f90 advec1D_BSL.f90 
	          advec2D_BSL.f90 
	          diag3D.f90 Pcross_section.f90 physics.f90)

SET(GYSELA_SRCS4 read_write.f90 output_saving.f90)

IF (NOT HDF5_FOUND)
   SET(GYSELA_SRCS ${GYSELA_SRCS1} ${GYSELA_SRCS3} ${GYSELA_SRCS4})
   ADD_DEFINITIONS(-DNOHDF5)
ELSE()
   SET(GYSELA_SRCS ${GYSELA_SRCS1} ${GYSELA_SRCS2} ${GYSELA_SRCS3} ${GYSELA_SRCS4})
ENDIF() 
ADD_LIBRARY(gysela_public STATIC ${GYSELA_SRCS})

ADD_DEPENDENCIES(gysela_public selalib)

ADD_EXECUTABLE(gysela diagnostics.f90 vlasov_solving.f90 gysela.f90)
TARGET_LINK_LIBRARIES(gysela gysela_public selalib ${HDF5_LIBRARIES})

ENABLE_TESTING()

INCLUDE(POEConfig)
SET(PROCS 16)
SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/../wk/data_test.nml)

ADD_MPI_TEST(test_gysela gysela ${PROCS} ${ARGS})

SET_TESTS_PROPERTIES(test_gysela PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
