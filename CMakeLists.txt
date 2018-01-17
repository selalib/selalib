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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(SeLaLib)

if (POLICY CMP0051)
  # CMake 3.1 and higher include generator expressions of the form
  # $<TARGETLIB:obj> in the SOURCES property.  These need to be
  # stripped everywhere that access the SOURCES property, so we just
  # defer to the OLD behavior of not including generator expressions
  # in the output for now.
  cmake_policy(SET CMP0051 OLD)
endif()
if (POLICY CMP0045)
  cmake_policy(SET CMP0045 OLD)
endif()

ENABLE_LANGUAGE(Fortran)

IF (${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_BINARY_DIR})
  MESSAGE(SEND_ERROR)
  MESSAGE("-- ERROR: in-tree-build not allowed.")
  MESSAGE("-- TRY: rm -f CMakeCache.txt; mkdir -p build; cd build; cmake ..; make")
  RETURN()
ENDIF()

OPTION(BUILD_HDF5             "Download and build the HDF5 library"    OFF)
OPTION(CLAPP                  "Use SELALIB within the CLAPP framework" OFF )
OPTION(OPENMP_ENABLED         "Enable use of OpenMP library"           OFF)
OPTION(MPI_ENABLED            "Enable use of MPI library"              ON )
OPTION(HDF5_ENABLED           "Use HDF5 library for data output"       ON )
OPTION(HDF5_PARALLEL_ENABLED  "Use Parallel HDF5"                      OFF)
OPTION(FFTW_ENABLED           "Enable use of FFTW library"             ON )
OPTION(PASTIX_ENABLED         "Enable PasTiX solver "                  OFF)
OPTION(MUMPS_ENABLED          "Enable MUMPS solver "                   OFF)
OPTION(SUITESPARSE_ENABLED    "Enable UMFpack solver"                  OFF)
OPTION(BUILD_SIMULATIONS      "Build selalib simulations"              ON )
OPTION(BUILD_GSL              "Build GSL and FGSL libraries"           OFF)
OPTION(BUILD_PARALUTION       "Build paralution library"               OFF)
OPTION(BUILD_SPRNG            "Build SPRNG libraries"                  OFF)
OPTION(USE_MKL                "Using Intel Math Kernel Library"        OFF)
OPTION(BUILD_THIRDPARTY       "Build third party simulations"          OFF)

INCLUDE(CMakeDependentOption)
CMAKE_DEPENDENT_OPTION(BUILD_MUDPACK "Build mudpack library" ON "BUILD_THIRDPARTY" ON)

SET(PROCS                 4   CACHE STRING "Number of MPI processes for tests")

IF(BUILD_THIRDPARTY)
    SET(HDF5_PARALLEL_ENABLED ON)
ENDIF(BUILD_THIRDPARTY)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  GET_FILENAME_COMPONENT(INSTALL_DIR "${CMAKE_SOURCE_DIR}/../usr" ABSOLUTE)
  SET(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH
    "Default installation directory" FORCE)
ELSE()
  SET(INSTALL_DIR ${CMAKE_INSTALL_PREFIX} )
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

INCLUDE(GitVersionString)

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

FILE(GLOB_RECURSE SLL_MACROS ${PROJECT_SOURCE_DIR}/src/*.h)
FILE(COPY ${SLL_MACROS} DESTINATION ${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/modules)
INCLUDE(Forcheck)

IF(MPI_ENABLED)
  INCLUDE(MPIConfig)
  INCLUDE(ProcessorCount)
  INCLUDE(POEConfig)
  INCLUDE(SLURMConfig)
ENDIF()

INCLUDE(CompilerConfigString)
INCLUDE(FortranConfig)
INCLUDE(HDF5Config)
INCLUDE(FFTConfig)
INCLUDE(DoxygenConfig)

# if BLASLAPACK_DIR was specified
IF(DEFINED BLASLAPACK_DIR)
  SET(BLASLAPACK_DIR ${BLASLAPACK_DIR} CACHE STRING "blas-lapack directory.")
  FIND_PACKAGE(BlasLapack)
ELSE()
  INCLUDE(LapackConfig)
ENDIF()

INCLUDE(PreprocessorTarget)
INCLUDE(CLAPP_framework)

INCLUDE(CTest)
ENABLE_TESTING()

MACRO(SLL_TEST _NAME)
  ADD_EXECUTABLE( test_${_NAME} "test_${_NAME}.F90")
  TARGET_LINK_LIBRARIES(test_${_NAME} ${ARGN})
  ADD_TEST(NAME ${_NAME} COMMAND test_${_NAME})
  SET_TESTS_PROPERTIES(${_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDMACRO(SLL_TEST)

ADD_SUBDIRECTORY(external)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(thirdparty)

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
  SET( INSTALL_DIR_SIMS ${INSTALL_DIR} )
  SET(CMAKE_INSTALL_PREFIX_SIMS ${INSTALL_DIR_SIMS} CACHE PATH
    "simulations install dir")
  ADD_SUBDIRECTORY( simulations )
ENDIF( BUILD_SIMULATIONS )

# CMake target for running C preprocessor on all source files used in build
ADD_PREPROCESSOR_TARGET()
ADD_FORCHECK_TARGET()

MARK_AS_ADVANCED(
  CMAKE_OSX_ARCHITECTURES
  CMAKE_OSX_DEPLOYMENT_TARGET
  CMAKE_OSX_SYSROOT
  CTEST_GIT_COMMAND
  DL_LIBRARY
  LSB_RELEASE_COMMAND
  OMPI_INFO_PROG
  PROCESSOR_COUNT
  ZLIB_LIBRARIES
  cmd_sys_pro PROCS SZ_LIBRARY 
  UNUSED_DUMMY_WARNING_ENABLED
  UNUSED_FUNCTION_WARNING_ENABLED
  FFTW_ROOT HDF5_ROOT
)
