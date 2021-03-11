add_executable(test_io test_io.F90)
target_link_libraries(test_io sll_file_io)
add_test(NAME file_io COMMAND test_io)
set_tests_properties(file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_ascii_io test_ascii_io.F90)
target_link_libraries(test_ascii_io sll_file_io)
add_test(NAME ascii_file_io COMMAND test_ascii_io)
set_tests_properties(ascii_file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_gnuplot test_gnuplot.F90)
target_link_libraries(test_gnuplot sll_file_io)
add_test(NAME gnuplot_file_io COMMAND test_gnuplot)
set_tests_properties(gnuplot_file_io PROPERTIES PASS_REGULAR_EXPRESSION
                                                "PASSED")

if(HDF5_ENABLED)
  add_executable(test_hdf5_io_serial test_hdf5_io_serial.F90)
  target_link_libraries(test_hdf5_io_serial sll_file_io)
  add_test(NAME hdf5_io_serial COMMAND test_hdf5_io_serial)
  set_tests_properties(hdf5_io_serial PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED" TIMEOUT 1)
endif(HDF5_ENABLED)

add_executable(test_io_tri_mesh test_io_tri_mesh.F90 m_generate_tri_mesh.F90)
target_link_libraries(test_io_tri_mesh sll_file_io)
