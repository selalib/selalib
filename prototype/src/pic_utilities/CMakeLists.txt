
add_library( sll_pic_utilities STATIC
             sll_pic_utilities.F90
             sll_representation_conversion.F90 
             sll_charge_to_density.F90
             sll_pic_particle_sort.F90 
             sll_quadtree.F90
             )

add_dependencies( sll_pic_utilities
                  sll_pic_accumulators
		  sll_assertion                  
		  sll_logical_meshes
                  sll_memory
		  sll_working_precision )

 
#add_executable(test_picutilities
#               unit_test.F90 ) 
#target_link_libraries( test_picutilities
#                       sll_pic_utilities
#		       sll_logical_meshes
#                       sll_assertion
#		       sll_memory
#		       sll_working_precision )
#
add_executable( test_particle_sort_2d
		unit_test_particle_sort.F90 )
target_link_libraries( test_particle_sort_2d 
                       sll_random_deviate_generators
		       sll_pic_utilities
		       sll_pic_particles
		       sll_logical_meshes      )

add_executable( test_quadtree
		    sll_unit_test_quadtree.F90 )

target_link_libraries(test_quadtree
			  sll_pic_utilities)
		    