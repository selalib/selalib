ADD_LIBRARY ( sll_hermite_interpolation STATIC
		sll_hermite_interpolation_2d.F90
			     	)

TARGET_LINK_LIBRARIES(sll_hermite_interpolation sll_utilities
                           sll_boundary_condition_descriptors
)

ADD_LIBRARY ( hex_hermite STATIC
              sll_hex_interpolation_hermite.F90
	      )

TARGET_LINK_LIBRARIES( hex_hermite sll_hermite_interpolation sll_hex_meshes )


# ADD_LIBRARY ( hex_poisson STATIC hex_poisson.F90)

# TARGET_LINK_LIBRARIES( hex_poisson sll_hex_meshes sll_hermite_interpolation)

# ADD_LIBRARY ( euler_2d_hex STATIC euler_2d_hex.F90)

# TARGET_LINK_LIBRARIES( euler_2d_hex )
# ADD_LIBRARY ( pivotbande STATIC pivotbande.F90)

# TARGET_LINK_LIBRARIES( pivotbande )


# first executable

ADD_EXECUTABLE(test_hermite_interpolation unit_test.F90)

TARGET_LINK_LIBRARIES(test_hermite_interpolation sll_hermite_interpolation
			                                sll_constants)

# second executable

ADD_EXECUTABLE(test_hex_hermite test_hex_hermite.F90)

TARGET_LINK_LIBRARIES(test_hex_hermite 
  sll_hermite_interpolation 
  hex_hermite 
  sll_hex_meshes 
  sll_constants)

# third executable

ADD_EXECUTABLE(test_hctc test_hctc.F90)

TARGET_LINK_LIBRARIES(test_hctc sll_constants sll_memory)

# fourth executable

ADD_EXECUTABLE(test_dimitrov test_dimitrov.F90)

TARGET_LINK_LIBRARIES(test_dimitrov sll_constants sll_memory)

# fifth executable

ADD_EXECUTABLE( test_hex_center_guide 
                test_hex_center_guide.F90 
                euler_2d_hex.F90
		../poisson_solvers/pivotbande.F90 )

TARGET_LINK_LIBRARIES( test_hex_center_guide 
                       sll_hermite_interpolation 
		       hex_hermite 
		       sll_hex_meshes 
		       hex_poisson
		       sll_constants )
		    
