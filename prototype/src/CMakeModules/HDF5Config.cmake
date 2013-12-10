SET(HDF5_ENABLED          ON              CACHE BOOL "Use HDF5 format for data output ")
SET(HDF5_PARALLEL_ENABLED OFF             CACHE BOOL "Use Parallel HDF5")
SET(HDF5_ROOT             $ENV{HDF5_ROOT} CACHE PATH "HDF5 location")

IF(NOT HDF5_FOUND AND HDF5_ENABLED)

   SET(HDF5_PATHS $ENV{HDF5_HOME}
                  ${HDF5_ROOT} 
                  $ENV{HDF5ROOT} 
                  /usr 
                  /usr/lib64/mpich2 
                  /usr/lib64/openmpi 
                  /usr/local 
                  /opt/local)

   FIND_PATH(HDF5_INCLUDE_DIR NAMES H5pubconf.h
   HINTS ${HDF5_PATHS} $ENV{HDF5_INCLUDEDIR} /usr/include/openmpi-x86_64 /usr/include/mpich2-x86_64 
   PATH_SUFFIXES / include hdf5/include 
   DOC "PATH to H5pubconf.h")

   FIND_PATH(HDF5_INCLUDE_DIR_FORTRAN NAMES hdf5.mod
   HINTS ${HDF5_PATHS} $ENV{HDF5_INCLUDEDIR} /usr/include/openmpi-x86_64 /usr/include/mpich2-x86_64 
   PATH_SUFFIXES / include hdf5/include include/fortran
   DOC "PATH to hdf5.mod")

   FIND_LIBRARY(HDF5_C_LIBRARY NAMES libhdf5.a hdf5
   HINTS ${HDF5_PATHS} $ENV{HDF5_LIBRARYDIR}
   PATH_SUFFIXES lib hdf5/lib lib/x86_64-linux-gnu
   DOC "PATH TO libhdf5")

   FIND_LIBRARY(HDF5_FORTRAN_LIBRARY NAMES libhdf5_fortran.a hdf5_fortran
   HINTS ${HDF5_PATHS} $ENV{HDF5_LIBRARYDIR}
   PATH_SUFFIXES lib hdf5/lib lib/x86_64-linux-gnu
   DOC "PATH TO libhdf5_fortran")

   FIND_LIBRARY(ZLIB_LIBRARIES NAMES z sz
                HINTS ${HDF5_PATHS} 
	          PATH_SUFFIXES lib hdf5/lib
	          DOC "PATH TO zip library")

   IF(HDF5_INCLUDE_DIR AND HDF5_INCLUDE_DIR_FORTRAN AND 
      HDF5_FORTRAN_LIBRARY AND ZLIB_LIBRARIES)
      SET(HDF5_FOUND YES)
   ENDIF()

   SET(HDF5_LIBRARIES ${HDF5_FORTRAN_LIBRARY} ${HDF5_C_LIBRARY} ${ZLIB_LIBRARIES})

   MESSAGE(STATUS "HDF5_INCLUDE_DIR:${HDF5_INCLUDE_DIR}")
   MESSAGE(STATUS "HDF5_INCLUDE_DIR_FORTRAN:${HDF5_INCLUDE_DIR_FORTRAN}")
   MESSAGE(STATUS "HDF5_LIBRARIES:${HDF5_LIBRARIES}")
   MESSAGE(STATUS "ZLIB_LIBRARIES:${ZLIB_LIBRARIES}")

ENDIF()


IF(HDF5_FOUND)

   MESSAGE(STATUS "HDF5 FOUND")

   set( HDF5_IS_PARALLEL FALSE )
   
   if( EXISTS "${HDF5_INCLUDE_DIR}/H5pubconf.h" )
      file( STRINGS "${HDF5_INCLUDE_DIR}/H5pubconf.h" 
          HDF5_HAVE_PARALLEL_DEFINE
          REGEX "HAVE_PARALLEL 1" )
      if( HDF5_HAVE_PARALLEL_DEFINE )
         set( HDF5_IS_PARALLEL TRUE )
      endif()
   endif()
   set( HDF5_IS_PARALLEL ${HDF5_IS_PARALLEL} CACHE BOOL
       "HDF5 library compiled with parallel IO support" )
   mark_as_advanced( HDF5_IS_PARALLEL )

   IF(HDF5_IS_PARALLEL) 
      MESSAGE(STATUS "HDF5 parallel supported")
   ELSE(HDF5_IS_PARALLEL)
      MESSAGE(STATUS "HDF5 parallel not supported")
   ENDIF()

   FIND_LIBRARY(GPFS_LIBRARY NAMES gpfs)
   IF(GPFS_LIBRARY)
      SET(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${GPFS_LIBRARY})
   ENDIF()

   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR_FORTRAN})
   MESSAGE(STATUS "HDF5_LIBRARIES:${HDF5_LIBRARIES}")

ELSE()

   MESSAGE(STATUS "Build SeLaLib without HDF5... binary output only for serial applications ")
   ADD_DEFINITIONS(-DNOHDF5)
   SET(HDF5_ENABLED OFF CACHE BOOL " " FORCE)
   SET(HDF5_LIBRARIES "")

ENDIF()
