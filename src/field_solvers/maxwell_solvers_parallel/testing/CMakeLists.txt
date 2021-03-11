add_executable(test_maxwell_2d_pstd_parallel test_maxwell_2d_pstd_parallel.F90)

target_link_libraries(test_maxwell_2d_pstd_parallel sll_maxwell_parallel
                      sll_file_io_parallel)
