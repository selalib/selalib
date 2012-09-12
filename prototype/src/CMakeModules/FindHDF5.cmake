#-----------------------------------------------------------------------------
# User Options
#-----------------------------------------------------------------------------
SET (HDF5_ENABLE_PARALLEL @HDF5_ENABLE_PARALLEL@)
SET (HDF5_BUILD_FORTRAN   @HDF5_BUILD_FORTRAN@)
SET (HDF5_ENABLE_F2003    @HDF5_ENABLE_F2003@)
SET (HDF5_BUILD_HL_LIB    @HDF5_BUILD_HL_LIB@)

find_path(HDF5_INCLUDE_DIRS NAMES hdf5.h hdf5.mod
	HINTS ${HDF5_ROOT}
	PATH_SUFFIXES include hdf5/include include/fortran
	DOC "PATH TO hdf5.h and hdf5.mod")

#SET(HDF5_INCLUDE_DIRS "@HDF5_INCLUDES_BUILD_TIME@" )

IF (HDF5_BUILD_FORTRAN)
  SET (HDF5_INCLUDE_DIR_FORTRAN "@CMAKE_Fortran_MODULE_DIRECTORY@" )
ENDIF (HDF5_BUILD_FORTRAN)

find_library(HDF5_HDF5_LIBRARY NAMES hdf5
	HINTS ${HDF5_ROOT}
	PATH_SUFFIXES lib hdf5/lib
	DOC "PATH TO libhdf5.dylib")

find_library(HDF5_HDF5_FORTRAN_LIBRARY NAMES hdf5_fortran
	HINTS ${HDF5_ROOT}
	PATH_SUFFIXES lib hdf5/lib
	DOC "PATH TO libhdf5_fortran.a")

find_library(HDF5_Z_LIBRARY NAMES z
	HINTS ${HDF5_ROOT}
	PATH_SUFFIXES lib hdf5/lib
	DOC "PATH TO libz.dylib")

set(HDF5_LIBRARIES @HDF5_HDF5_FORTRAN_LIBRARY@;@HDF5_HDF5_LIBRARY@;@HDF5_Z_LIBRARY@)

IF ( HDF5_INCLUDE_DIRS         AND
     HDF5_HDF5_LIBRARY         AND
     HDF5_HDF5_FORTRAN_LIBRARY AND
     HDF5_Z_LIBRARY                )
  set(HDF5_FOUND YES)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
ENDIF()

IF(HDF5_ENABLE_PARALLEL) 
   MESSAGE(STATUS "HDF5 parallel supported")
ELSE(HDF5_ENABLE_PARALLEL)
   MESSAGE(STATUS "HDF5 parallel not supported")
ENDIF()
IF(HDF5_BUILD_FORTRAN)   
   MESSAGE(STATUS "HDF5 was compiled with fortran on")
ENDIF()
IF(HDF5_BUILD_HL_LIB)    
   MESSAGE(STATUS "HDF5 was compiled with high level on")
ENDIF()
IF(HDF5_ENABLE_F2003)
   MESSAGE (STATUS "HDF5 FORTRAN 2003 Standard enabled")
ENDIF()
