add_executable( test_pic_visu_parallel test_pic_visu_parallel.F90 )
target_link_libraries( test_pic_visu_parallel sll_pic_visu_parallel)
ADD_MPI_TEST( pic_visu_parallel test_pic_visu_parallel 4 "" )
