ADD_EXECUTABLE( test_mesh_calculus_2d test_mesh_calculus_2d.F90 )
TARGET_LINK_LIBRARIES( test_mesh_calculus_2d
    sll_mesh_calculus sll_file_io sll_memory ${HDF5_LIBRARIES} )
