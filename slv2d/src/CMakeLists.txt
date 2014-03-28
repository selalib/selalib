SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release." )

IF(OPTIONS_FILE)
INCLUDE(${OPTIONS_FILE})
ENDIF()

PROJECT(SLV2D)
ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


MESSAGE( STATUS "SLV2D_BINARY_DIR: ${SLV2D_BINARY_DIR}" )
MESSAGE( STATUS "SLV2D_SOURCE_DIR: ${SLV2D_SOURCE_DIR}" )
IF (CMAKE_PROJECT_NAME MATCHES "SLV2D")
   SET(SLL_SOURCE_DIR "${SLV2D_SOURCE_DIR}/../../prototype/src")
ENDIF()
MESSAGE( STATUS "SLL_SOURCE_DIR: ${SLL_SOURCE_DIR}" )

SET (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake 
                       ${SLL_SOURCE_DIR}/CMakeModules
                       ${CMAKE_MODULE_PATH} )

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
INCLUDE_DIRECTORIES(${SLL_SOURCE_DIR}/package)

SET(ENABLE_PERF OFF CACHE BOOL "Enable perf library linking")
IF(ENABLE_PERF)
  FIND_LIBRARY(OOPERF_LIBRARY ooperf /home/rop/lib/PERF)
  FIND_LIBRARY(PFM_LIBRARY pfm /home/rop/lib/PERF)
  SET(PERF_LIBRARIES ${OOPERF_LIBRARY} ${PFM_LIBRARY})
  LINK_LIBRARIES("${PERF_LIBRARIES} -lstdc++")
ENDIF(ENABLE_PERF)

MESSAGE(STATUS "${CMAKE_PROJECT_NAME}:CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
INCLUDE(F90FLAGS)
MESSAGE(STATUS "${CMAKE_PROJECT_NAME}:CMAKE_Fortran_FLAGS:${CMAKE_Fortran_FLAGS}")

IF(NOT LAPACK_FOUND)
   INCLUDE(${SLL_SOURCE_DIR}/CMakeModules/LapackConfig.cmake)
ENDIF(NOT LAPACK_FOUND)

SET(SILO_ENABLED OFF CACHE BOOL "Use SILO Library from LLNL ")
SET(FFTW_ENABLED ON CACHE BOOL "Use fftw library " )
SET(MPI_MODULE_ENABLED ON CACHE BOOL "Use MPI library " )
SET(HDF5_ENABLED ON CACHE BOOL "Enable HDF5 ")
SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "Enable HDF5 parallel")

INCLUDE(MPIConfig)
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})


INCLUDE(HDF5Config)
IF(HDF5_IS_PARALLEL) 
   MESSAGE(STATUS "HDF5 parallel supported")
ELSE()
    MESSAGE(SEND_ERROR "HDF5 parallel not supported")
    MESSAGE(SEND_ERROR "./configure --enable-fortran --enable-parallel CC=mpicc FC=mpif90")
ENDIF()

IF(SILO_ENABLED) 
   FIND_PACKAGE(SILO  REQUIRED)
   IF(SILO_FOUND)
      ADD_DEFINITIONS(-D_SILO)
      INCLUDE_DIRECTORIES(${SILO_INCLUDE_DIRS})
      SET(CMAKE_EXE_LINKER_FLAGS "-lsupc++")
   ENDIF()
ENDIF()

IF(FFTW_ENABLED AND NOT FFTW_FOUND)
   FIND_PACKAGE(FFTW)
   IF(FFTW_FOUND)
      ADD_DEFINITIONS(-D_FFTW)
      INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIRS})
   ENDIF()
ENDIF()



SET(FFT_LIB "FFTW")
SET(FFTPACK_LIBRARIES "dfftpack")

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)

IF (CMAKE_PROJECT_NAME MATCHES "SLV2D")
   ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/../../prototype/src/fftpack 
                    ${CMAKE_BINARY_DIR}/dfftpack)
   ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/../../prototype/src/package 
                    ${CMAKE_BINARY_DIR}/package)
ENDIF()

SET(SLV2D_SRC silo_module.F90
              diagnostiques.F90
              geometry.F90
              maxwell2dfdtd.F90
              fft_module.F90
              poisson2dpp_seq.F90
              poisson2d_periodic.F90
              precision.F90
              splinenn.F90
              splinepp.F90
              splinepx.F90
              splinepy.F90
              transpose.F90
              vlasov1d.F90
              vlasov2d.F90 
              sll_vlasov4d.F90
              vlasov4d_plot.F90)

ADD_LIBRARY(slv2d STATIC ${SLV2D_SRC})
ADD_DEPENDENCIES(slv2d fftpack selalib)

MESSAGE(STATUS "${CMAKE_PROJECT_NAME}:SILO_LIBRARIES:${SILO_LIBRARIES}")
MESSAGE(STATUS "${CMAKE_PROJECT_NAME}:FFTPACK_LIBRARIES:${FFTPACK_LIBRARIES}")
MESSAGE(STATUS "${CMAKE_PROJECT_NAME}:FFTW_LIBRARIES:${FFTW_LIBRARIES}")
MESSAGE(STATUS "${CMAKE_PROJECT_NAME}:HDF5_LIBRARIES:${HDF5_LIBRARIES}")

SET(SLV2D_CSL2D_SRC silo_module.F90
              diagnostiques.F90
              geometry.F90
              maxwell2dfdtd.F90
              poisson2dpp_seq.F90
              fft_module.F90
              poisson2d_periodic.F90
              precision.F90
              splinenn.F90
              splinepp.F90
              csl2dpp.F90
              splinepx.F90
              splinepy.F90
              transpose.F90
              vlasov1d.F90
              vlasov2d_csl2d.F90 
              sll_vlasov4d.F90
              vlasov4d_plot.F90)

ADD_LIBRARY(slv2d_csl2d STATIC ${SLV2D_CSL2D_SRC})

ADD_DEPENDENCIES(slv2d_csl2d selalib-mpi dfftpack)


SET(SLV2D_NON_UNIF_V_SRC silo_module.F90
              diagnostiques.F90
              geometry.F90
              maxwell2dfdtd.F90
              poisson2dpp_seq.F90
              fft_module.F90
              poisson2d_periodic.F90
              precision.F90
              splinenn.F90
              splinepp.F90
              interp_non_unif_pp.F90
              splinepx.F90
              splinepy.F90
              transpose.F90
              vlasov1d.F90
              vlasov2d_non_unif_v.F90 
              sll_vlasov4d.F90
              vlasov4d_plot.F90)

ADD_LIBRARY(slv2d_non_unif_v STATIC ${SLV2D_NON_UNIF_V_SRC})

ADD_DEPENDENCIES(slv2d_non_unif_v selalib-mpi dfftpack)


SET(SLV2D_DK_SRC silo_module.F90
              diagnostiques.F90
              geometry.F90
              maxwell2dfdtd.F90
              poisson2dpp_seq.F90
              fft_module.F90
              poisson2d_periodic.F90
              precision.F90
              splinenn.F90
              splinepp.F90
              splinepx.F90
              splinepy.F90
              transpose.F90
              vlasov1d.F90
              vlasov2d_dk.F90 
              sll_vlasov4d.F90
              vlasov4d_plot.F90)

ADD_LIBRARY(slv2d_dk STATIC ${SLV2D_DK_SRC})

ADD_DEPENDENCIES(slv2d_dk selalib-mpi dfftpack)





MESSAGE(STATUS "SLV2D:SILO_LIBRARIES:${SILO_LIBRARIES}")
MESSAGE(STATUS "SLV2D:FFTPACK_LIBRARIES:${FFTPACK_LIBRARIES}")
MESSAGE(STATUS "SLV2D:FFTW_LIBRARIES:${FFTW_LIBRARIES}")
MESSAGE(STATUS "SLV2D:HDF5_LIBRARIES:${HDF5_LIBRARIES}")

SET(LIBS selalib 
         ${SILO_LIBRARIES} 
         ${FFTPACK_LIBRARIES} 
         ${FFTW_LIBRARIES} 
         ${HDF5_LIBRARIES})

ADD_EXECUTABLE(vp4d_transpose vp4d_transpose.F90)
TARGET_LINK_LIBRARIES(vp4d_transpose slv2d ${LIBS} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(vp4d_csl2d vp4d_csl2d.F90)
TARGET_LINK_LIBRARIES(vp4d_csl2d slv2d_csl2d ${LIBS} ${LAPACK_LIBRARIES})


ADD_EXECUTABLE(vp4d_dk vp4d_dk.F90)
TARGET_LINK_LIBRARIES(vp4d_dk slv2d_dk ${LIBS} ${LAPACK_LIBRARIES})


ADD_EXECUTABLE(vp2d_keen_non_unif_v vp2d_keen_non_unif_v.F90)
TARGET_LINK_LIBRARIES(vp2d_keen_non_unif_v slv2d_non_unif_v ${LIBS} ${LAPACK_LIBRARIES})


ADD_EXECUTABLE(vm4d_transpose vm4d_transpose.F90)
TARGET_LINK_LIBRARIES(vm4d_transpose slv2d ${LIBS} ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(vp4d_transpose_cubic vp4d_cubic.F90)
TARGET_LINK_LIBRARIES(vp4d_transpose_cubic slv2d ${LIBS})

ADD_EXECUTABLE(vp2d_keen vp2d_keen.F90)
TARGET_LINK_LIBRARIES(vp2d_keen slv2d ${LIBS} ${LAPACK_LIBRARIES})

ENABLE_TESTING()

INCLUDE(POEConfig)
SET(PROCS 4)
FILE(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../input/slv2d.nml" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

SET(ARGS "slv2d.nml")

ADD_MPI_TEST(test_vp4d       vp4d_transpose       ${PROCS} ${ARGS})
ADD_MPI_TEST(test_vm4d       vm4d_transpose       ${PROCS} ${ARGS})
ADD_MPI_TEST(test_vp4d_cubic vp4d_transpose_cubic ${PROCS} ${ARGS})

SET_TESTS_PROPERTIES(test_vp4d          PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(test_vm4d          PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
SET_TESTS_PROPERTIES(test_vp4d_cubic    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

