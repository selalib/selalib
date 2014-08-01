ADD_LIBRARY (sll_hermite_interpolation STATIC
		sll_hermite_interpolation_2d.F90)

TARGET_LINK_LIBRARIES(sll_hermite_interpolation sll_utilities  
                           hex_mesh
                           sll_boundary_condition_descriptors)

# first executable

ADD_EXECUTABLE(test_hermite_interpolation unit_test.F90)

TARGET_LINK_LIBRARIES(test_hermite_interpolation sll_hermite_interpolation
			                                sll_constants)

# second executable

ADD_EXECUTABLE(test_hex_hermite test_hex_hermite.F90)

TARGET_LINK_LIBRARIES(test_hex_hermite sll_hermite_interpolation 
				sll_constants)
