SET(HDF5_ENABLED ON CACHE BOOL "Use HDF5 format for data output ")
SET(HDF5_PARALLEL_ENABLED OFF CACHE BOOL "Use Parallel HDF5")

IF(NOT HDF5_FOUND)

   SET(HDF5_PATHS $ENV{HDF5_HOME}
                  $ENV{HDF5_ROOT} 
                  /usr 
                  /usr/lib64/mpich2 
                  /usr/lib64/openmpi 
                  /usr/local 
                  /opt/local)

   FIND_PATH(HDF5_INCLUDE_DIRS NAMES hdf5.mod
   HINTS ${HDF5_PATHS} /usr/include/openmpi-x86_64 /usr/include/mpich2-x86_64
   PATH_SUFFIXES / include hdf5/include include/fortran
   DOC "PATH to hdf5.mod")

   FIND_LIBRARY(HDF5_C_LIBRARY NAMES hdf5
   HINTS ${HDF5_PATHS} 
   PATH_SUFFIXES lib hdf5/lib
   DOC "PATH TO libhdf5")

   FIND_LIBRARY(HDF5_FORTRAN_LIBRARY NAMES hdf5_fortran
   HINTS ${HDF5_PATHS} 
   PATH_SUFFIXES lib hdf5/lib
   DOC "PATH TO libhdf5_fortran")

   FIND_LIBRARY(ZLIB_LIBRARIES NAMES z
                HINTS ${HDF5_PATHS} 
	        PATH_SUFFIXES lib hdf5/lib
	        DOC "PATH TO libz.dylib")

   IF(HDF5_INCLUDE_DIRS AND HDF5_FORTRAN_LIBRARY AND ZLIB_LIBRARIES)
      SET(HDF5_FOUND YES)
   ENDIF()

   SET(HDF5_LIBRARIES ${HDF5_FORTRAN_LIBRARY} ${HDF5_C_LIBRARY} ${ZLIB_LIBRARIES})

   MESSAGE(STATUS "HDF5_INCLUDE_DIRS:${HDF5_INCLUDE_DIRS}")
   MESSAGE(STATUS "HDF5_LIBRARIES:${HDF5_LIBRARIES}")
   MESSAGE(STATUS "ZLIB_LIBRARIES:${ZLIB_LIBRARIES}")

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
   #INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS}/fortran)
   MESSAGE(STATUS "HDF5_LIBRARIES:${HDF5_LIBRARIES}")

ELSE()

   MESSAGE(STATUS "Build SeLaLib without HDF5... binary output only for serial applications ")
   ADD_DEFINITIONS(-DNOHDF5)
   SET(HDF5_ENABLED OFF CACHE BOOL " " FORCE)
   SET(HDF5_LIBRARIES "")

ENDIF()
