add_executable( test_pic_visu_parallel test_pic_visu_parallel.F90 )
target_link_libraries( test_pic_visu_parallel sll_pic_visu_parallel)
ADD_MPI_TEST( pic_visu_parallel test_pic_visu_parallel 4 "" )
add_test(pic_visu_parallel_check_files ${CMAKE_COMMAND} -E compare_files ${CMAKE_CURRENT_BINARY_DIR}/visu_pic_parallel_test_file_0000.xmf ${CMAKE_CURRENT_SOURCE_DIR}/visu_pic_parallel_test_file_0000.xmf)
