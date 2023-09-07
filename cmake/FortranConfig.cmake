# Determine how-to install the modules. CMAKE_BINARY_DIR is the directory in
# which the make command is invoked.
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Add the modules directory to the list of include directories
# INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

get_filename_component(Fortran_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)
message(STATUS "CMAKE_Fortran_COMPILER_ID:${CMAKE_Fortran_COMPILER_ID}")

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

  exec_program(
    ${CMAKE_Fortran_COMPILER} ARGS
    "-v"
    OUTPUT_VARIABLE source_path)
  message(STATUS "${source_path}")
  string(REGEX MATCH "1[0-9]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION
               ${source_path})
  set(CMAKE_Fortran_FLAGS_RELEASE "-nowarn -O3 -xHost -ip -fpic")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-g -O0 -check all,noarg_temp_created -fpe0 -traceback -ftrapuv -fpic")
  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-shared-intel")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES IBM)

  set(CMAKE_Fortran_FLAGS_DEBUG "-qextname=flush -qxlf2003=polymorphic")
  set(CMAKE_Fortran_FLAGS_RELEASE
      "-qnosave -qextname=flush -qxlf2003=polymorphic")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)

  exec_program(
    ${CMAKE_Fortran_COMPILER} ARGS
    "--version"
    OUTPUT_VARIABLE source_path)
  string(REGEX MATCH "([0-9]+)\\.[0-9]+\\.[0-9]+" Fortran_COMPILER_VERSION
               ${source_path})

  add_definitions(-DGFORTRAN)
  set(CMAKE_Fortran_FLAGS_RELEASE
      "-std=f2008 -ffree-line-length-none -fstack-arrays -O3 -fPIC  -w ")
  if(NOT APPLE)
    set(CMAKE_Fortran_FLAGS_RELEASE
        "${CMAKE_Fortran_FLAGS_RELEASE} -march=native")
  endif()
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-std=f2008 -ffree-line-length-none \
  -fstack-arrays -O0 -g -fbacktrace -Werror=intrinsics-std -Wall \
  -pedantic -Wconversion-extra -Wuninitialized \
  -fcheck=array-temps,bounds,do,pointer,recursion \
  -ffpe-trap=invalid,zero,overflow") # -Wno-integer-division -Werror")

  set(UNUSED_FUNCTION_WARNING_ENABLED
      OFF
      CACHE BOOL "Add -Wunused-function flag to gfortran")
  if(NOT UNUSED_FUNCTION_WARNING_ENABLED)
    set(CMAKE_Fortran_FLAGS_DEBUG
        "${CMAKE_Fortran_FLAGS_DEBUG} -Wno-unused-function")
  endif()

  set(UNUSED_DUMMY_WARNING_ENABLED
      OFF
      CACHE BOOL "Add -Wunused-dummy-argument flag to gfortran")
  if(NOT UNUSED_DUMMY_WARNING_ENABLED)
    set(CMAKE_Fortran_FLAGS_DEBUG
        "${CMAKE_Fortran_FLAGS_DEBUG} -Wno-unused-dummy-argument")
  endif()

  if(Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10")
    set(CMAKE_Fortran_FLAGS_DEBUG
        "${CMAKE_Fortran_FLAGS_DEBUG} -fallow-argument-mismatch")
    set(CMAKE_Fortran_FLAGS_RELEASE
        "${CMAKE_Fortran_FLAGS_RELEASE} -fallow-argument-mismatch")
  endif()

else()

  message(SEND_ERROR "NO KNOWN FORTRAN COMPILER FOUND")

endif()

message(
  STATUS "Fortran ${Fortran_COMPILER_NAME} version ${Fortran_COMPILER_VERSION}")

# --- enable fully user-defineable compiler flags
if(FORCE_Fortran_FLAGS_RELEASE)
  set(CMAKE_Fortran_FLAGS_RELEASE "${FORCE_Fortran_FLAGS_RELEASE}")
endif()

set(ADDITIONAL_COMPILER_FLAGS
    ""
    CACHE STRING "The user can define additional compiler flags here")
set(CMAKE_Fortran_FLAGS_DEBUG
    "${CMAKE_Fortran_FLAGS_DEBUG} ${ADDITIONAL_COMPILER_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE
    "${CMAKE_Fortran_FLAGS_RELEASE} ${ADDITIONAL_COMPILER_FLAGS}")

if(OPENMP_ENABLED)
  find_package(OpenMP_Fortran)
  set(CMAKE_Fortran_FLAGS_DEBUG
      "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_RELEASE
      "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
endif()

mark_as_advanced(CLEAR CMAKE_Fortran_COMPILER)
mark_as_advanced(CLEAR CMAKE_C_COMPILER)
mark_as_advanced(CLEAR CMAKE_CXX_COMPILER)
