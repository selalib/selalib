ADD_LIBRARY(sll_meshes STATIC sll_meshes_base.F90 sll_cartesian_meshes.F90)

TARGET_LINK_LIBRARIES( sll_meshes sll_working_precision sll_memory)

ADD_EXECUTABLE( test_meshes unit_test.F90)
TARGET_LINK_LIBRARIES( test_meshes sll_meshes)

ADD_TEST(NAME meshes            COMMAND test_meshes)
SET_TESTS_PROPERTIES(meshes PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_LIBRARY(sll_hex_meshes STATIC sll_meshes_base.F90 sll_hex_meshes.F90)

TARGET_LINK_LIBRARIES(sll_hex_meshes sll_file_io
                                     sll_utilities
                                     sll_working_precision 
		                     sll_memory
			 	     sll_constants)

ADD_EXECUTABLE( test_hex_meshes test_hex_meshes.F90)
TARGET_LINK_LIBRARIES( test_hex_meshes sll_hex_meshes)
