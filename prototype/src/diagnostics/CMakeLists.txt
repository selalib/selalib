# use all sll_*F90 files

add_custom_command(
	OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/sll_hdf5_io.F90
	COMMAND m4
	ARGS ${CMAKE_CURRENT_SOURCE_DIR}/sll_hdf5_io.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/sll_hdf5_io.F90
	VERBATIM
)

add_custom_command(
	OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/sll_ascii_io.F90
	COMMAND m4
	ARGS ${CMAKE_CURRENT_SOURCE_DIR}/sll_ascii_io.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/sll_ascii_io.F90
	VERBATIM
)

add_custom_command(
	OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/sll_binary_io.F90
	COMMAND m4
	ARGS ${CMAKE_CURRENT_SOURCE_DIR}/sll_binary_io.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/sll_binary_io.F90
	VERBATIM
)

file(GLOB DIAGNOS_SRC "sll_*.F90")
add_library(sll_diagnostics STATIC ${DIAGNOS_SRC}
	sll_hdf5_io.F90 sll_binary_io.F90 sll_ascii_io.F90
)
target_link_libraries( sll_diagnostics ${HDF5_LIBRARIES} )

add_executable(test_diagnostics unit_test.F90)
target_link_libraries(test_diagnostics sll_diagnostics sll_memory sll_assertion)
