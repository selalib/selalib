file(GLOB IO_SRC "sll_*.F90")
#set(IO_SRC sll_ascii_io.F90 sll_binary_io.F90 sll_gnuplot.F90 sll_hdf5_io.F90 sll_io.F90 sll_xdmf.F90 sll_xml_io.F90)
add_library(sll_low_level_file_io STATIC ${IO_SRC})

IF(HDF5_FOUND AND HDF5_ENABLED)
	target_link_libraries(sll_low_level_file_io ${HDF5_LIBRARIES})
ENDIF()

add_executable(test_io unit_test.F90)
target_link_libraries(test_io
	sll_low_level_file_io
	sll_memory
	sll_assertion
)
