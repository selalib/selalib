ADD_EXECUTABLE( sim2d_gc_splines_hex simulation_2d_guiding_center_hexagonal.F90)

TARGET_LINK_LIBRARIES( sim2d_gc_splines_hex 
		       pivotbande
		       euler_2d_hex
                       sll_box_splines
		       sll_hex_meshes 
		       hex_poisson
		       sll_constants )
		    
