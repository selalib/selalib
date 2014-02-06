ADD_LIBRARY( sll_scalar_field_multipatch STATIC 
             sll_scalar_field_2d_multipatch.F90 )		

TARGET_LINK_LIBRARIES( sll_scalar_field_multipatch
                       sll_coordinate_transformations_multipatch
		       sll_logical_meshes_multipatch )

ADD_EXECUTABLE( test_scalar_field_multipatch_2d 
                unit_test_2d.F90 )

TARGET_LINK_LIBRARIES( test_scalar_field_multipatch_2d
                       sll_scalar_field_multipatch)

