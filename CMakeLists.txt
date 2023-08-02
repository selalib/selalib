if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      ${CMAKE_BUILD_TYPE}
      CACHE STRING "Choose the type of build, options are: Debug Release.")
else()
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build, options are: Debug Release.")
endif()

if(OPTIONS_FILE)
  include(${OPTIONS_FILE})
  message("OPTIONS_FILE ARE READ")
endif()

cmake_minimum_required(VERSION 3.6.0)
project(SeLaLib)

enable_language(Fortran)

if(${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_BINARY_DIR})
  message(SEND_ERROR)
  message("-- ERROR: in-tree-build not allowed.")
  message(
    "-- TRY: rm -f CMakeCache.txt; mkdir -p build; cd build; cmake ..; make")
  return()
endif()

option(OPENMP_ENABLED "Enable use of OpenMP library" OFF)
option(MPI_ENABLED "Enable use of MPI library" ON)
option(HDF5_ENABLED "Use HDF5 library for data output" ON)
option(HDF5_PARALLEL_ENABLED "Use Parallel HDF5" ON)
option(FFTW_ENABLED "Enable use of FFTW library" ON)
option(PASTIX_ENABLED "Enable PasTiX solver " OFF)
option(MUMPS_ENABLED "Enable MUMPS solver " OFF)
option(SUITESPARSE_ENABLED "Enable UMFpack solver" OFF)
option(BUILD_SIMULATIONS "Build selalib simulations" ON)
option(BUILD_FGSL "Build GSL and FGSL libraries" OFF)
option(USE_MKL "Using Intel Math Kernel Library" OFF)

include(CMakeDependentOption)
cmake_dependent_option(BUILD_MUDPACK "Build mudpack library" ON
                       "BUILD_THIRDPARTY" ON)

set(PROCS
    4
    CACHE STRING "Number of MPI processes for tests")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  get_filename_component(INSTALL_DIR "${CMAKE_SOURCE_DIR}/../usr" ABSOLUTE)
  set(CMAKE_INSTALL_PREFIX
      ${INSTALL_DIR}
      CACHE PATH "Default installation directory" FORCE)
else()
  set(INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
endif()

find_package(PythonInterp 3)
set(PYTHON3_FOUND FALSE)
message(STATUS "PYTHON_EXECUTABLE:${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_VERSION_STRING:${PYTHON_VERSION_STRING}")
if(PYTHON_EXECUTABLE AND PYTHON_VERSION_STRING MATCHES "3.*")
  string(COMPARE GREATER ${PYTHON_VERSION_STRING} "3.0.0" PYTHON3_FOUND)
endif()
message(STATUS "PYTHON3_FOUND:${PYTHON3_FOUND}")
mark_as_advanced(CLEAR PYTHON_EXECUTABLE)

# Specify the output for binary
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add new directory to find the FIND<PACKAGE>.cmake files. CMAKE_SOURCE_DIR is
# the directory where the project sources are located. CMAKE_MODULE_PATH is the
# list of directories in which to search for CMake modules.
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DDEBUG)
endif()

include(GitVersionString)

# Define folder to be deleted by make clean
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                    ${CMAKE_BINARY_DIR}/modules/)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                    ${CMAKE_BINARY_DIR}/bin/)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                                    ${CMAKE_BINARY_DIR}/include/)

add_custom_target(
  clear
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
  COMMENT "Delete all .a, .mod, .o and folders /modules, /bin"
  VERBATIM)

file(GLOB_RECURSE SLL_MACROS ${PROJECT_SOURCE_DIR}/src/*.h)
file(COPY ${SLL_MACROS} DESTINATION ${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/modules)
include(Forcheck)

# INCLUDE(CompilerConfigString)
include(FortranConfig)

if(MPI_ENABLED)
  include(MPIConfig)
  if(MPI_FOUND)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  endif(MPI_FOUND)
  include(ProcessorCount)
  include(POEConfig)
  include(SLURMConfig)
endif()

include(HDF5Config)
include(ZFPConfig)
include(FFTConfig)
include(DoxygenConfig)

find_package(BLAS)
find_package(LAPACK)

include(PreprocessorTarget)
# Note: The following include must be placed after the compiler flags are fixed.
include(CompilerConfigString)

include(CTest)
enable_testing()

macro(SLL_TEST _NAME)
  add_executable(test_${_NAME} "test_${_NAME}.F90")
  target_link_libraries(test_${_NAME} ${ARGN})
  add_test(NAME ${_NAME} COMMAND test_${_NAME})
  set_tests_properties(${_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endmacro(SLL_TEST)

add_subdirectory(external)
add_subdirectory(src)

message(STATUS "############################################################")
message(STATUS "####")
message(STATUS "#### CMAKE_BUILD_TYPE:${CMAKE_BUILD_TYPE}")
message(STATUS "#### Fortran_COMPILER_NAME:${Fortran_COMPILER_NAME}")
message(STATUS "#### CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
message(STATUS "#### Fortran FLAGS ")
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  message(STATUS "#### ${CMAKE_Fortran_FLAGS_DEBUG}")
else()
  message(STATUS "#### ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()
message(STATUS "####")
message(STATUS "############################################################")

add_custom_target(
  clean_outputs
  COMMAND find . -name *.mtv -delete
  COMMAND find . -name *.gnu -delete
  COMMAND find . -name *.dat -delete
  COMMAND find . -name *.xmf -delete
  COMMAND find . -name *.h5 -delete
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Delete all .gnu, .dat, .mtv, .h5, .xmf "
  VERBATIM)

# Create single library 'libselalib.a' from smaller libraries in src
add_subdirectory(package)

# Simulations
if(BUILD_SIMULATIONS)
  set(INSTALL_DIR_SIMS ${INSTALL_DIR})
  set(CMAKE_INSTALL_PREFIX_SIMS
      ${INSTALL_DIR_SIMS}
      CACHE PATH "simulations install dir")
  add_subdirectory(simulations)
endif(BUILD_SIMULATIONS)

# CMake target for running C preprocessor on all source files used in build
# ADD_PREPROCESSOR_TARGET() ADD_FORCHECK_TARGET()

mark_as_advanced(
  CMAKE_OSX_ARCHITECTURES
  CMAKE_OSX_DEPLOYMENT_TARGET
  CMAKE_OSX_SYSROOT
  CTEST_GIT_COMMAND
  DL_LIBRARY
  LSB_RELEASE_COMMAND
  OMPI_INFO_PROG
  PROCESSOR_COUNT
  ZLIB_LIBRARIES
  cmd_sys_pro
  PROCS
  SZ_LIBRARY
  UNUSED_DUMMY_WARNING_ENABLED
  UNUSED_FUNCTION_WARNING_ENABLED)
