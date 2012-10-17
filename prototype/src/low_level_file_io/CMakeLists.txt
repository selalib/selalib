#FILE(GLOB IO_SRC "sll_*.F90")
SET(IO_SRC sll_ascii_io.F90 
           sll_binary_io.F90 
           sll_gnuplot.F90 
           sll_io.F90 
           sll_hdf5_io.F90
           sll_xml_io.F90
           sll_xdmf.F90 )


ADD_LIBRARY(sll_low_level_file_io STATIC ${IO_SRC})

ADD_DEPENDENCIES(sll_low_level_file_io sll_assertion sll_working_precision sll_utilities)

ADD_EXECUTABLE(test_io unit_test.F90)

TARGET_LINK_LIBRARIES(sll_low_level_file_io ${HDF5_LIBRARIES})

TARGET_LINK_LIBRARIES(test_io sll_low_level_file_io
                              sll_constants
                              sll_memory
                              sll_assertion
                              sll_working_precision)

IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)

   ADD_DEFINITIONS(-DHDF5_PARALLEL)
   SET(IO_SRC ${IO_SRC} sll_hdf5_io_parallel.F90
                        sll_xdmf_parallel.F90 )
   ADD_LIBRARY(sll_low_level_file_io_parallel STATIC ${IO_SRC})
   ADD_DEPENDENCIES(sll_low_level_file_io_parallel sll_utilities 
                                                   sll_collective)
   TARGET_LINK_LIBRARIES(sll_low_level_file_io_parallel ${HDF5_LIBRARIES} )
   ADD_EXECUTABLE(test_io_parallel unit_test_parallel.F90)
   TARGET_LINK_LIBRARIES(test_io_parallel sll_low_level_file_io_parallel
                                          sll_remap
                                          sll_collective
	                                  sll_utilities
	                                  sll_constants
	                                  sll_memory
	                                  sll_assertion
	                                  sll_working_precision)
ENDIF()
