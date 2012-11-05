
# add the cache entry HDF5_ENABLED for enable/disable hdf5
set(HDF5_ENABLED ON CACHE BOOL "Use HDF5 format for data output ")
set(HDF5_PARALLEL_ENABLED OFF CACHE BOOL "Use Parallel HDF5")
# search for HDF5
IF(HDF5_ENABLED)
	IF( DEFINED ENV{HDF5_ROOT} )
		set(HDF5_ROOT $ENV{HDF5_ROOT})
	ENDIF()
        #IF( DEFINED ENV{HDF5_HOME} )
	#	set(HDF5_ROOT $ENV{HDF5_HOME})
	#	message(STATUS "HDF5_HOME FOUND")
	#ENDIF()
	find_package(HDF5 QUIET)
	IF(HDF5_FOUND)
		message(STATUS "HDF5 FOUND")
	ELSE()
		add_definitions(-DNOHDF5)
		set(HDF5_ENABLED OFF CACHE BOOL " " FORCE)
		set(HDF5_LIBRARIES "")
	ENDIF()
ELSE(HDF5_ENABLED)
	add_definitions(-DNOHDF5)
ENDIF(HDF5_ENABLED)
##########################################################
