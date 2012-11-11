# add the cache entry HDF5_ENABLED for enable/disable hdf5
SET(HDF5_ENABLED ON CACHE BOOL "Use HDF5 format for data output ")
SET(HDF5_PARALLEL_ENABLED OFF CACHE BOOL "Use Parallel HDF5")

IF($ENV{HOSTNAME} MATCHES "hpc-f0*")
   SET(HDF5_USE_STATIC_LIBRARIES YES)
   SET(HDF5_ROOT "/home/math/navaro/local/bin")
ELSE()
   SET(HDF5_ROOT $ENV{HDF5_ROOT})
   FIND_PACKAGE(HDF5 REQUIRED Fortran)
ENDIF()


IF(NOT HDF5_FOUND)

   MESSAGE(STATUS "CMake did not find your HDF5 installation")

   FIND_PATH(HDF5_INCLUDE_DIRS NAMES hdf5.h
   HINTS ${HDF5_ROOT}/../include /usr/include /usr/lib64/mpich2/include /usr/lib64/openmpi/include /usr/local/include
   PATH_SUFFIXES include hdf5/include
   DOC "PATH TO hdf5.h")

   FIND_PATH(HDF5_INCLUDE_DIR_FORTRAN NAMES hdf5.mod
   HINTS ${HDF5_ROOT}/../include /usr/include /usr/lib64/mpich2/include /usr/lib64/openmpi/include /usr/local/include
   PATH_SUFFIXES include hdf5/include include/fortran
   DOC "PATH to hdf5.mod")

   FIND_LIBRARY(HDF5_HDF5_LIBRARY NAMES hdf5
   HINTS ${HDF5_ROOT}/../lib /usr/lib /usr/lib64/mpich2/lib /usr/lib64/openmpi/lib /usr/local/lib
   PATH_SUFFIXES lib hdf5/lib
   DOC "PATH TO libhdf5")

   FIND_LIBRARY(HDF5_HDF5_FORTRAN_LIBRARY NAMES hdf5_fortran
   HINTS ${HDF5_ROOT}/../lib /usr/lib /usr/lib64/mpich2/lib /usr/lib64/openmpi/lib /usr/local/lib
   PATH_SUFFIXES lib hdf5/lib
   DOC "PATH TO libhdf5_fortran")

   FIND_PACKAGE(ZLIB)

   IF (ZLIB_FOUND)
      SET(HDF5_LIBRARIES @HDF5_HDF5_FORTRAN_LIBRARY@;@HDF5_HDF5_LIBRARY@ ${ZLIB_LIBRARIES})
      SET(HDF5_Z_LIBRARY ${ZLIB_LIBRARIES})
   ELSE()
      FIND_LIBRARY(HDF5_Z_LIBRARY NAMES z
	           HINTS ${HDF5_ROOT}
	           PATH_SUFFIXES lib hdf5/lib
	           DOC "PATH TO libz.dylib")
   ENDIF()

   IF ( HDF5_INCLUDE_DIRS         AND
        HDF5_HDF5_LIBRARY         AND
        HDF5_HDF5_FORTRAN_LIBRARY AND
        HDF5_Z_LIBRARY )

     MESSAGE(STATUS "Ok we have everything we need to link with HDF5")
     
     SET(HDF5_FOUND YES)
     SET(HDF5_LIBRARIES @HDF5_HDF5_FORTRAN_LIBRARY@;@HDF5_HDF5_LIBRARY@;@HDF5_Z_LIBRARY@)

   ENDIF()

ENDIF()


IF(HDF5_FOUND)

   MESSAGE(STATUS "HDF5 FOUND")
   SET(HDF5_ENABLE_PARALLEL @HDF5_ENABLE_PARALLEL@)
   SET(HDF5_BUILD_FORTRAN   @HDF5_BUILD_FORTRAN@)
   SET(HDF5_ENABLE_F2003    @HDF5_ENABLE_F2003@)
   SET(HDF5_BUILD_HL_LIB    @HDF5_BUILD_HL_LIB@)
   
   IF(HDF5_ENABLE_PARALLEL) 
      MESSAGE(STATUS "HDF5 parallel supported")
   ELSE(HDF5_ENABLE_PARALLEL)
      MESSAGE(STATUS "HDF5 parallel not supported")
   ENDIF()
   IF(HDF5_BUILD_FORTRAN)   
      MESSAGE(STATUS "HDF5 was compiled with fortran on")
   ELSE(HDF5_BUILD_FORTRAN)   
      MESSAGE(STATUS "HDF5 was compiled with fortran off")
   ENDIF()
   IF(HDF5_BUILD_HL_LIB)    
      MESSAGE(STATUS "HDF5 was compiled with high level on")
   ELSE(HDF5_BUILD_HL_LIB)    
      MESSAGE(STATUS "HDF5 was compiled with high level off")
   ENDIF()
   IF(HDF5_ENABLE_F2003)
      MESSAGE (STATUS "HDF5 FORTRAN 2003 Standard enabled")
   ELSE(HDF5_ENABLE_F2003)
      MESSAGE (STATUS "HDF5 FORTRAN 2003 Standard disabled")
   ENDIF()

   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
   INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS}/fortran)  

ELSE()

   MESSAGE(STATUS "Build SeLaLib without HDF5... binary output only for serial applications ")
   ADD_DEFINITIONS(-DNOHDF5)
   SET(HDF5_ENABLED OFF CACHE BOOL " " FORCE)
   SET(HDF5_LIBRARIES "")

ENDIF()
