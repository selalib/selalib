SET(CMAKE_BUILD_TYPE Release CACHE STRING 
    "Choose the type of build, options are: Debug Release." )
PROJECT(vlasov4d)
ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(SLL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../prototype/src)

SET(CMAKE_MODULE_PATH ${SLL_SOURCE_DIR}/CMakeModules 
                       ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake
                       ${CMAKE_MODULE_PATH})

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")


MESSAGE(STATUS "Fortran COMPILER:${CMAKE_Fortran_COMPILER}")
INCLUDE(F90FLAGS)
MESSAGE(STATUS "Fortran FLAGS:${CMAKE_Fortran_FLAGS_DEBUG}")

INCLUDE(${SLL_SOURCE_DIR}/CMakeModules/LapackConfig.cmake)
FIND_PACKAGE(FFTW)
FIND_PACKAGE(MPI)

SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})


SET(FFTW_ENABLED ON CACHE BOOL "Enable FFTW" FORCE)
SET(HDF5_ENABLED ON CACHE BOOL "Enable HDF5" FORCE)
SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "Enable HDF5 parallel" FORCE)

INCLUDE(HDF5Config)

IF(HDF5_IS_PARALLEL) 
   MESSAGE(STATUS "HDF5 parallel supported")
ELSE()
   MESSAGE(SEND_ERROR "HDF5 parallel not supported")
ENDIF()

SET(FFT_LIB "FFTW")
ADD_DEFINITIONS(-DFFTW)

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)

IF (CMAKE_PROJECT_NAME MATCHES "vlasov4d")
  INCLUDE_DIRECTORIES(${SLL_SOURCE_DIR}/package)
  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/../../prototype/src/package 
                 ${CMAKE_BINARY_DIR}/package)
  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/../../prototype/src/fftpack 
                 ${CMAKE_BINARY_DIR}/fftpack)
ENDIF()

SET(LIBS selalib ${FFTW_LIBRARIES} ${HDF5_LIBRARIES})

ADD_LIBRARY(selalib_vlasov4d STATIC sll_vlasov4d_spectral.F90
                                    sll_vlasov4d_spectral_charge.F90
                                    sll_vlasov4d_maxwell.F90
                                    sll_vlasov4d_poisson.F90
                                    sll_vlasov4d_polar.F90
                                    sll_vlasov4d_base.F90)

ADD_DEPENDENCIES(selalib_vlasov4d selalib)

ADD_EXECUTABLE(vp4d          vp4d_remapper.F90)
ADD_EXECUTABLE(vm4d          vm4d_remapper.F90)
ADD_EXECUTABLE(vm4d_spectral vm4d_spectral.F90)
ADD_EXECUTABLE(vp4d_polar    vp4d_polar.F90)
ADD_EXECUTABLE(vm4d_spectral_charge vm4d_spectral_charge.F90 init_functions.F90)

TARGET_LINK_LIBRARIES(vp4d          selalib_vlasov4d ${LIBS})
TARGET_LINK_LIBRARIES(vm4d          selalib_vlasov4d ${LIBS})
TARGET_LINK_LIBRARIES(vm4d_spectral selalib_vlasov4d ${LIBS})
TARGET_LINK_LIBRARIES(vp4d_polar    selalib_vlasov4d ${LIBS})
TARGET_LINK_LIBRARIES(vm4d_spectral_charge selalib_vlasov4d ${LIBS})

ADD_SUBDIRECTORY(vm4d-csl)
