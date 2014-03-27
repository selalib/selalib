ADD_LIBRARY(sll_logical_meshes STATIC 
  sll_logical_meshes.F90)

TARGET_LINK_LIBRARIES( sll_logical_meshes sll_working_precision sll_memory)


ADD_LIBRARY(hex_mesh STATIC hex_mesh.F90)

TARGET_LINK_LIBRARIES( hex_mesh sll_working_precision sll_memory)

ADD_EXECUTABLE( test_logical_meshes unit_test.F90)

TARGET_LINK_LIBRARIES( test_logical_meshes sll_logical_meshes)

ADD_EXECUTABLE( test_hex_mesh test_hex_mesh.F90)

TARGET_LINK_LIBRARIES( test_hex_mesh hex_mesh)
