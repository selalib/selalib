ADD_EXECUTABLE( sim2d_gc_hex_splines simulation_2d_guiding_center_hexagonal_splines.F90)

TARGET_LINK_LIBRARIES( sim2d_gc_hex_splines
		       pivotbande
		       euler_2d_hex
                       sll_box_splines
		       sll_hex_meshes 
		       hex_poisson
		       sll_constants )
		    

ADD_EXECUTABLE( sim2d_gc_hex_hermite simulation_2d_guiding_center_hexagonal_hermite.F90)

TARGET_LINK_LIBRARIES( sim2d_gc_hex_hermite
		       pivotbande
		       euler_2d_hex
                       sll_box_splines
		       sll_hex_meshes 
		       hex_poisson
		       sll_constants )
