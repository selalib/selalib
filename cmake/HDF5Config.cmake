set(HDF5_ROOT
    $ENV{HDF5_ROOT}
    CACHE PATH "HDF5 location")

set( HDF5_USE_STATIC_LIBRARIES ON )
set( HDF5_PREFER_PARALLEL TRUE)
set( HDF5_FIND_DEBUG TRUE )

find_package( HDF5 COMPONENTS Fortran )

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

if(HDF5_FOUND)

  macro(CHECK_HDF5_DEPS HDF5_HAVE_STRING HDF5_HAVE_BOOL)
    file(STRINGS "${HDF5_INCLUDE_DIR}/H5pubconf.h" HDF5_HAVE_DEFINE
         REGEX ${HDF5_HAVE_STRING})
    if(HDF5_HAVE_DEFINE)
      set(${HDF5_HAVE_BOOL} TRUE)
    else()
      set(${HDF5_HAVE_BOOL} FALSE)
    endif()
  endmacro(CHECK_HDF5_DEPS)

  if(EXISTS "${HDF5_INCLUDE_DIR}/H5pubconf.h")
    check_hdf5_deps("HAVE_LIBPTHREAD 1" HDF5_HAVE_LIBPTHREAD)
    check_hdf5_deps("HAVE_GPFS 1" HDF5_HAVE_GPFS)
    check_hdf5_deps("HAVE_LIBDL 1" HDF5_HAVE_LIBDL)
    check_hdf5_deps("HAVE_LIBSZ 1" HDF5_HAVE_LIBSZ)
  endif()

  set(HDF5_HAVE_LIBPTHREAD
      ${HDF5_HAVE_LIBPTHREAD}
      CACHE BOOL "HDF5 library compiled with pthread library")
  mark_as_advanced(HDF5_HAVE_LIBPTHREAD)
  if(HDF5_HAVE_LIBPTHREAD)
    find_library(PTHREAD_LIBRARY NAMES pthread)
    set(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${PTHREAD_LIBRARY})
  endif()

  set(HDF5_HAVE_GPFS
      ${HDF5_HAVE_GPFS}
      CACHE BOOL "HDF5 library compiled with GPFS")
  mark_as_advanced(HDF5_HAVE_GPFS)
  if(HDF5_HAVE_GPFS)
    find_library(GPFS_LIBRARY NAMES gpfs)
    set(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${GPFS_LIBRARY})
  endif()

  set(HDF5_HAVE_LIBDL
      ${HDF5_HAVE_LIBDL}
      CACHE BOOL "HDF5 library compiled with LIBDL")
  mark_as_advanced(HDF5_HAVE_LIBDL)
  if(HDF5_HAVE_LIBDL)
    find_library(DL_LIBRARY NAMES dl)
    set(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${DL_LIBRARY})
  endif()

  set(HDF5_HAVE_LIBSZ
      ${HDF5_HAVE_LIBSZ}
      CACHE BOOL "HDF5 library compiled with LIBSZ")
  mark_as_advanced(HDF5_HAVE_LIBSZ)
  if(HDF5_HAVE_LIBSZ)
    find_library(SZ_LIBRARY NAMES sz)
    set(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${SZ_LIBRARY})
  endif()

else()

  message( STATUS "Build SeLaLib without HDF5... binary output only")
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
