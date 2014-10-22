ADD_LIBRARY(sll_meshes_base STATIC sll_meshes_base.F90)

TARGET_LINK_LIBRARIES( sll_meshes_base sll_working_precision sll_memory)

ADD_LIBRARY(sll_logical_meshes STATIC sll_logical_meshes.F90)

TARGET_LINK_LIBRARIES( sll_logical_meshes sll_meshes_base
                                          sll_working_precision 
					  sll_memory)

ADD_LIBRARY(sll_hex_meshes STATIC sll_hex_meshes.F90)

TARGET_LINK_LIBRARIES(sll_hex_meshes sll_meshes_base
                                     sll_file_io
                                     sll_utilities
                                     sll_working_precision 
		                     sll_memory
			 	     sll_constants)

ADD_EXECUTABLE( test_logical_meshes unit_test.F90)

TARGET_LINK_LIBRARIES( test_logical_meshes sll_logical_meshes)

ADD_EXECUTABLE( test_hex_meshes test_hex_meshes.F90)

TARGET_LINK_LIBRARIES( test_hex_meshes sll_hex_meshes)
