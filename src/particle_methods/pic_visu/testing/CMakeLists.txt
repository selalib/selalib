add_executable(test_pic_visu test_pic_visu.F90 biot_savart.F90)
target_link_libraries(test_pic_visu sll_pic_visu)
set(PROCS 1)
set(ARGS " ")
add_mpi_test(pic_visu test_pic_visu ${PROCS} ${ARGS})
set(pic_visu PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_pic_viewer test_pic_viewer.F90 biot_savart.F90)
target_link_libraries(test_pic_viewer sll_pic_visu sll_poisson_solvers
                      sll_file_io sll_coordinate_transformations)
set(PROCS 1)
set(ARGS " ")
add_mpi_test(pic_viewer test_pic_viewer ${PROCS} ${ARGS})
set(pic_visu PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
