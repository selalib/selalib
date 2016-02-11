IF(DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} 
    CACHE STRING "Choose the type of build, options are: Debug Release.")
ELSE()
  SET(CMAKE_BUILD_TYPE Debug 
    CACHE STRING "Choose the type of build, options are: Debug Release.")
ENDIF()

IF(OPTIONS_FILE)
  INCLUDE(${OPTIONS_FILE})
  MESSAGE("OPTIONS_FILE ARE READ")
ENDIF()

PROJECT(SeLaLib)

ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

SET(OPENMP_ENABLED        OFF CACHE BOOL   "Enable use of OpenMP library")
SET(MPI_ENABLED           ON  CACHE BOOL   "Enable use of MPI library")
SET(HDF5_ENABLED          ON  CACHE BOOL   "Use HDF5 library for data output")
SET(HDF5_PARALLEL_ENABLED OFF CACHE BOOL   "Use Parallel HDF5")
SET(FFTW_ENABLED          ON  CACHE BOOL   "Enable use of FFTW library")
SET(PASTIX_ENABLED        OFF CACHE BOOL   "Build PasTiX solvers ")
SET(SUITESPARSE_ENABLED   OFF CACHE BOOL   "Enable UMFpack solvers")
SET(BUILD_SIMULATIONS     ON  CACHE BOOL   "Build selalib simulations")
SET(BUILD_GSL             OFF CACHE BOOL   "Build GSL and FGSL libraries")
SET(BUILD_PARALUTION      OFF CACHE BOOL   "Build paralution library")
SET(BUILD_SPRNG           OFF CACHE BOOL   "Build SPRNG libraries")
SET(PROCS                 4   CACHE STRING "Number of MPI processes for tests")
SET(USE_MKL               OFF CACHE BOOL   "Using Intel Math Kernel Library")
SET(BUILD_MUDPACK         OFF CACHE BOOL   "Build mudpack library")

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  GET_FILENAME_COMPONENT(INSTALL_DIR "${CMAKE_SOURCE_DIR}/usr" ABSOLUTE)
  SET(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH
    "Default installation directory" FORCE)
ENDIF()

FIND_PACKAGE(PythonInterp 3)
SET(PYTHON3_FOUND FALSE)
MESSAGE(STATUS "PYTHON_EXECUTABLE:${PYTHON_EXECUTABLE}")
MESSAGE(STATUS "PYTHON_VERSION_STRING:${PYTHON_VERSION_STRING}")
IF(PYTHON_EXECUTABLE AND PYTHON_VERSION_STRING MATCHES "3.*")
  STRING(COMPARE GREATER ${PYTHON_VERSION_STRING} "3.0.0" PYTHON3_FOUND)
ENDIF()
MESSAGE(STATUS "PYTHON3_FOUND:${PYTHON3_FOUND}")
MARK_AS_ADVANCED(CLEAR PYTHON_EXECUTABLE)

# Specify the output for binary
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add new directory to find the FIND<PACKAGE>.cmake files. CMAKE_SOURCE_DIR is
# the directory where the project sources are located. CMAKE_MODULE_PATH is 
# the list of directories in which to search for CMake modules.
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )

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
  COMMAND find . -name *.xmf -delete
  COMMAND find . -name *.h5 -delete
  COMMAND find . -name *.gnu -delete
  COMMAND find . -name *.dat -delete
  COMMAND rm -rf bin
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Delete all .a, .mod, .o and folders /modules, /bin" VERBATIM)

INCLUDE(FortranConfig)
INCLUDE(FFTConfig)
INCLUDE(LapackConfig)
INCLUDE(DoxygenConfig)
INCLUDE(PreprocessorTarget)

IF(MPI_ENABLED)
  INCLUDE(MPIConfig)
  INCLUDE(ProcessorCount)
  INCLUDE(POEConfig)
ENDIF()

INCLUDE(CTest)
INCLUDE(HDF5Config)

IF(PASTIX_ENABLED)
  FIND_PACKAGE(SCOTCH)
  FIND_PACKAGE(PASTIX)
ELSE()
  IF(SUITESPARSE_ENABLED)
    FIND_PACKAGE(SUITESPARSE)
  ENDIF(SUITESPARSE_ENABLED)
ENDIF(PASTIX_ENABLED)

MACRO(SLL_TEST _NAME)
  ADD_EXECUTABLE( test_${_NAME} "test_${_NAME}.F90")
  TARGET_LINK_LIBRARIES(test_${_NAME} ${ARGN})
  ADD_TEST(NAME ${_NAME} COMMAND test_${_NAME})
  SET(${_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDMACRO(SLL_TEST)

ENABLE_TESTING()

FILE(GLOB_RECURSE SLL_MACROS ${PROJECT_SOURCE_DIR}/src/*.h)
FILE(COPY ${SLL_MACROS} DESTINATION ${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/modules)

ADD_SUBDIRECTORY(external)
ADD_SUBDIRECTORY(src)

IF(BUILD_PARALUTION)
  INCLUDE(ExternalProject)
  SET(PARALUTION_VERSION 1.1.0 CACHE STRING "Paralution version number")
  EXTERNALPROJECT_ADD( paralution
     URL  http://www.paralution.com/downloads/paralution-${PARALUTION_VERSION}.tar.gz
     CMAKE_ARGS -DSUPPORT_CUDA=OFF -DBUILD_EXAMPLES=OFF
     INSTALL_COMMAND ""
  )
ENDIF()

IF($ENV{PFUNIT})
  INCLUDE(ExternalProject)
  EXTERNALPROJECT_ADD(
    pFUnit
    PREFIX pFunit
    GIT_REPOSITORY git://git.code.sf.net/p/pfunit/code
    CMAKE_ARGS -DINSTALL_PATH=${CMAKE_BINARY_DIR})
ENDIF()

IF(BUILD_SPRNG)
  INCLUDE(BuildSPRNG)
ENDIF()


MESSAGE(STATUS "############################################################")
MESSAGE(STATUS "####")
MESSAGE(STATUS "#### CMAKE_BUILD_TYPE:${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "#### Fortran_COMPILER_NAME:${Fortran_COMPILER_NAME}")
MESSAGE(STATUS "#### CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
MESSAGE(STATUS "#### Fortran FLAGS ")
IF(CMAKE_BUILD_TYPE MATCHES "Debug") 
   MESSAGE(STATUS "#### ${CMAKE_Fortran_FLAGS_DEBUG}")
ELSE()
   MESSAGE(STATUS "#### ${CMAKE_Fortran_FLAGS_RELEASE}")
ENDIF()
MESSAGE(STATUS "####")
MESSAGE(STATUS "############################################################")

ADD_CUSTOM_TARGET( clean_outputs
   COMMAND find . -name *.mtv -delete
   COMMAND find . -name *.gnu -delete
   COMMAND find . -name *.dat -delete
   COMMAND find . -name *.xmf -delete
   COMMAND find . -name *.h5  -delete
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   COMMENT "Delete all .gnu, .dat, .mtv, .h5, .xmf " VERBATIM)

# Create single library 'libselalib.a' from smaller libraries in src
ADD_SUBDIRECTORY( package )

# Simulations
IF( BUILD_SIMULATIONS )
  GET_FILENAME_COMPONENT( INSTALL_DIR_SIMS "${CMAKE_SOURCE_DIR}/usr" ABSOLUTE )
  SET(CMAKE_INSTALL_PREFIX_SIMS ${INSTALL_DIR_SIMS} CACHE PATH 
    "simulations install dir")
  ADD_SUBDIRECTORY( simulations )
ENDIF( BUILD_SIMULATIONS )

# CMake target for running C preprocessor on all source files used in build
ADD_PREPROCESSOR_TARGET()
