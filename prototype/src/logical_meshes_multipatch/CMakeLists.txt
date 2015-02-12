ADD_LIBRARY(sll_logical_meshes_multipatch STATIC 
  sll_logical_meshes_multipatch.F90)

TARGET_LINK_LIBRARIES( sll_logical_meshes_multipatch
 sll_logical_meshes
 sll_memory
 sll_assertion
 sll_working_precision )

ADD_EXECUTABLE( test_logical_meshes_multipatch unit_test.F90)

TARGET_LINK_LIBRARIES( test_logical_meshes_multipatch 
  sll_logical_meshes_multipatch
  sll_assertion )
