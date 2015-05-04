ADD_LIBRARY(sll_meshes_multipatch STATIC 
  sll_meshes_multipatch.F90)

TARGET_LINK_LIBRARIES( sll_meshes_multipatch
 sll_meshes
 sll_memory
 sll_assert
 sll_working_precision )


# Ctest
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_meshes_multipatch unit_test.F90)
  TARGET_LINK_LIBRARIES( test_meshes_multipatch sll_meshes_multipatch sll_assert )
  ADD_TEST( NAME meshes_multipatch COMMAND test_meshes_multipatch )
  SET_TESTS_PROPERTIES( meshes_multipatch PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)
