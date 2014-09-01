ADD_LIBRARY ( sll_hermite_interpolation STATIC
		sll_hermite_interpolation_2d.F90
			     	)

TARGET_LINK_LIBRARIES(sll_hermite_interpolation sll_utilities
                           sll_boundary_condition_descriptors
)

ADD_LIBRARY ( hex_hermite STATIC
                hex_interpolation_hermite.F90
			     	)

TARGET_LINK_LIBRARIES( hex_hermite sll_hermite_interpolation)

#TARGET_LINK_LIBRARIES( hex_mesh )



# first executable

ADD_EXECUTABLE(test_hermite_interpolation unit_test.F90)

TARGET_LINK_LIBRARIES(test_hermite_interpolation sll_hermite_interpolation
			                                sll_constants)

# second executable

ADD_EXECUTABLE(test_hex_hermite test_hex_hermite.F90)

TARGET_LINK_LIBRARIES(test_hex_hermite 
  sll_hermite_interpolation 
  hex_hermite 
  hex_mesh 
  sll_constants)

# third executable

ADD_EXECUTABLE(test_hctc test_hctc.F90)

TARGET_LINK_LIBRARIES(test_hctc sll_constants)


# fourth executable

ADD_EXECUTABLE(test_dimitrov test_dimitrov.F90)

TARGET_LINK_LIBRARIES(test_dimitrov sll_constants)