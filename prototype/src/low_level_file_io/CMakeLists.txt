#FILE(GLOB IO_SRC "sll_*.F90")
SET(IO_SRC sll_ascii_io.F90 
           sll_binary_io.F90 
           sll_gnuplot.F90 
           sll_io.F90 
           sll_xdmf.F90 
           sll_xml_io.F90)


IF(HDF5_FOUND AND HDF5_ENABLED)
   IF(HDF5_ENABLE_PARALLEL)
      SET(IO_SRC ${IO_SRC} sll_hdf5_io_parallel.F90)
   ELSE(HDF5_ENABLE_PARALLEL)
      SET(IO_SRC ${IO_SRC} sll_hdf5_io.F90)
   ENDIF(HDF5_ENABLE_PARALLEL)
ELSE()
   SET(HDF5_LIBRARIES "")
ENDIF()

ADD_LIBRARY(sll_low_level_file_io STATIC ${IO_SRC})
ADD_DEPENDENCIES(sll_low_level_file_io sll_assertion sll_working_precision)
TARGET_LINK_LIBRARIES(sll_low_level_file_io ${HDF5_LIBRARIES} )

ADD_EXECUTABLE(test_io unit_test.F90)
TARGET_LINK_LIBRARIES(test_io sll_low_level_file_io
	                      sll_constants
	                      sll_memory
	                      sll_assertion
	                      sll_working_precision)
