add_executable(test_mesh_calculus_2d test_mesh_calculus_2d.F90)
target_link_libraries(test_mesh_calculus_2d sll_mesh_calculus sll_file_io
                      sll_memory ${HDF5_LIBRARIES})
