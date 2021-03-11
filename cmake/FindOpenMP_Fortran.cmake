# This module can be used to detect OpenMP support in a compiler. If the
# compiler supports OpenMP, the flags required to compile with openmp support
# are set.
#
# This module was modified from the standard FindOpenMP module to find Fortran
# flags.
#
# The following variables are set: OpenMP_Fortran_FLAGS - flags to add to the
# Fortran compiler for OpenMP support.  In general, you must use these at both
# compile- and link-time. OMP_NUM_PROCS - the max number of processors available
# to OpenMP

# =============================================================================
# Copyright 2009 Kitware, Inc. Copyright 2008-2009 André Rigland Brodtkorb
# &lt;Andre.Brodtkorb@ifi.uio.no&gt;
#
# Distributed under the OSI-approved BSD License (the "License"); see
# accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# License for more information.
# =============================================================================
# (To distribute this file outside of CMake, substitute the full License text
# for the above reference.)
#
# [YG] 27.01.2017: Change Intel flag from -openmp to -qopenmp

include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

# Corner case: If both the environment variables are defined the build may fail
# during OpenMP detection.
if(DEFINED ENV{OMP_STACKSIZE} AND DEFINED ENV{KMP_STACKSIZE})
  unset(ENV{KMP_STACKSIZE})
endif()

set(OpenMP_Fortran_FLAG_CANDIDATES
    # Intel
    "-qopenmp"
    # Gnu
    "-fopenmp"
    # Empty, if compiler automatically accepts openmp
    " "
    # Sun
    "-xopenmp"
    # HP
    "+Oopenmp"
    # IBM XL C/c++
    "-qsmp"
    # Portland Group
    "-mp"
    # Microsoft Visual Studio
    "/openmp"
    # Intel windows
    "/Qopenmp")

if(DEFINED OpenMP_Fortran_FLAGS)
  set(OpenMP_Fortran_FLAG_CANDIDATES)
endif(DEFINED OpenMP_Fortran_FLAGS)

# check fortran compiler. also determine number of processors
foreach(FLAG ${OpenMP_Fortran_FLAG_CANDIDATES})
  set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_FLAGS "${FLAG}")
  unset(OpenMP_FLAG_DETECTED CACHE)
  message(STATUS "Try OpenMP Fortran flag = [${FLAG}]")
  file(
    WRITE
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMP.f90"
    "
program TestOpenMP
 use omp_lib
 write(*,'(I2)',ADVANCE='NO') omp_get_num_procs()
end program TestOpenMP
")
  set(MACRO_CHECK_FUNCTION_DEFINITIONS
      "-DOpenMP_FLAG_DETECTED ${CMAKE_REQUIRED_FLAGS}")
  try_run(
    OpenMP_RUN_FAILED OpenMP_FLAG_DETECTED ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranOpenMP.f90
    COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
    CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
    COMPILE_OUTPUT_VARIABLE OUTPUT
    RUN_OUTPUT_VARIABLE OMP_NUM_PROCS_INTERNAL)
  if(OpenMP_FLAG_DETECTED)
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
         "Determining if the Fortran compiler supports OpenMP passed with "
         "the following output:\n${OUTPUT}\n\n")
    set(OpenMP_FLAG_DETECTED 1)
    if(OpenMP_RUN_FAILED)
      message(FATAL_ERROR "OpenMP found, but test code did not run")
    endif(OpenMP_RUN_FAILED)
    set(OMP_NUM_PROCS
        ${OMP_NUM_PROCS_INTERNAL}
        CACHE STRING "Number of processors OpenMP may use" FORCE)
    set(OpenMP_Fortran_FLAGS_INTERNAL "${FLAG}")
    break()
  else()
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
         "Determining if the Fortran compiler supports OpenMP failed with "
         "the following output:\n${OUTPUT}\n\n")
    set(OpenMP_FLAG_DETECTED 0)
  endif(OpenMP_FLAG_DETECTED)
endforeach(FLAG ${OpenMP_Fortran_FLAG_CANDIDATES})

set(OpenMP_Fortran_FLAGS
    "${OpenMP_Fortran_FLAGS_INTERNAL}"
    CACHE STRING "Fortran compiler flags for OpenMP parallization")

# handle the standard arguments for FIND_PACKAGE
find_package_handle_standard_args(OpenMP_Fortran DEFAULT_MSG
                                  OpenMP_Fortran_FLAGS)

mark_as_advanced(OpenMP_Fortran_FLAGS)
