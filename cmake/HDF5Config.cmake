set(HDF5_ROOT
    $ENV{HDF5_ROOT}
    CACHE PATH "HDF5 location")

set( HDF5_USE_STATIC_LIBRARIES ON )
set( HDF5_PREFER_PARALLEL TRUE)

find_package( HDF5 COMPONENTS C Fortran REQUIRED )

message(STATUS "HDF5 FOUND : ${HDF5_FOUND}")
message(STATUS "HDF5 VERSION : ${HDF5_VERSION}")

if(HDF5_IS_PARALLEL)
   message(STATUS "HDF5 parallel supported")
   add_definitions(-DHDF5_PARALLEL)
else(HDF5_IS_PARALLEL)
   message(STATUS "HDF5 parallel not supported")
   message(STATUS "try to find the parallel version somewhere else")
endif()

if(NOT HDF5_IS_PARALLEL)

  if(MPI_Fortran_LIBRARIES)
    list(GET MPI_Fortran_LIBRARIES 0 HEAD)
    get_filename_component(MPI_LIB_ROOT ${HEAD} PATH)
  else()
    set(MPI_LIB_ROOT "")
  endif()

  set(HDF5_PATHS
      ${HDF5_ROOT}
      ${MPI_Fortran_INCLUDE_PATH}
      $ENV{HDF5_HOME}
      $ENV{HDF5ROOT}
      $ENV{HDF5_BASE}
      $ENV{HDF5_ROOT_DIR}
      $ENV{HDF5_DIR}
      $ENV{SZIP_LIB}
      $ENV{SZIP_LIBDIR}
      $ENV{HDF5_BASE}
      /usr
      /usr/include/hdf5/openmpi
      /usr/include/hdf5/mpich
      /usr/local
      /opt/local)

  find_path(
    HDF5_INCLUDE_DIR
    NAMES H5pubconf.h
    HINTS ${HDF5_PATHS} $ENV{HDF5_INCLUDEDIR} $ENV{HDF5_INC_DIR}
    PATH_SUFFIXES / include
    DOC "PATH to H5pubconf.h")

  find_path(
    HDF5_INCLUDE_DIR_FORTRAN
    NAMES hdf5.mod
    HINTS ${HDF5_PATHS} $ENV{HDF5_INCLUDEDIR} $ENV{HDF5_INC_DIR}
    PATH_SUFFIXES / include hdf5/include include/fortran
    DOC "PATH to hdf5.mod")

  find_library(
    HDF5_C_LIBRARY
    NAMES libhdf5.a hdf5_openmpi hdf5_mpich hdf5
    HINTS ${MPI_LIB_ROOT} ${HDF5_PATHS} $ENV{HDF5_LIBRARYDIR} $ENV{HDF5_LIB_DIR}
    PATH_SUFFIXES lib hdf5/lib lib/x86_64-linux-gnu
    DOC "PATH TO libhdf5")

  find_library(
    HDF5_FORTRAN_LIBRARY
    NAMES libhdf5_fortran.a hdf5_openmpi_fortran hdf5_mpich_fortran hdf5_fortran
    HINTS ${MPI_LIB_ROOT} ${HDF5_PATHS} $ENV{HDF5_LIBRARYDIR} $ENV{HDF5_LIB_DIR}
    PATH_SUFFIXES lib hdf5/lib lib/x86_64-linux-gnu
    DOC "PATH TO libhdf5_fortran")

  find_library(
    ZLIB_LIBRARIES
    NAMES z sz
    HINTS ${HDF5_PATHS}
    PATH_SUFFIXES lib hdf5/lib ENV${SZIP_LIB}
    DOC "PATH TO zip library")

  set(HDF5_LIBRARIES ${HDF5_FORTRAN_LIBRARY} ${HDF5_C_LIBRARY}
                     ${ZLIB_LIBRARIES})

  set(HDF5_Fortran_INCLUDE_DIRS  ${HDF5_INCLUDE_DIR_FORTRAN})
  set(HDF5_Fortran_LIBRARIES  ${HDF5_LIBRARIES})

  if(DEFINED HDF5_FORTRAN_LIBRARY)
    set(HDF5_FOUND YES)
  endif()

endif()

if(NOT HDF5_FOUND)

  message(
    STATUS
      "Build SeLaLib without HDF5... binary output only for serial applications "
  )
  add_definitions(-DNOHDF5)
  set(HDF5_ENABLED
      OFF
      CACHE BOOL " " FORCE)
  set(HDF5_LIBRARIES "")

endif()

if(HDF5_ENABLED
   AND HDF5_IS_PARALLEL
   AND NOT MPI_ENABLED)

  message(STATUS "HDF5 is PARALLEL and needs MPI, please set MPI_ENABLED")
  message(STATUS "HDF5 is set to OFF")
  set(HDF5_ENABLED
      OFF
      CACHE BOOL " " FORCE)
  add_definitions(-DNOHDF5)

endif()

message(STATUS "##########################################################")
message(STATUS "HDF5_INCLUDE_DIRS:${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5_Fortran_INCLUDE_DIRS:${HDF5_Fortran_INCLUDE_DIRS}")
message(STATUS "HDF5_Fortran_LIBRARIES:${HDF5_Fortran_LIBRARIES}")
message(STATUS "##########################################################")
