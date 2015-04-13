
add_library( sll_pic_utilities STATIC
             sll_pic_utilities.F90
             sll_representation_conversion.F90 
             sll_charge_to_density.F90
             sll_pic_particle_sort.F90
             sll_particle_method_descriptors.F90)

target_link_libraries( sll_pic_utilities
  sll_pic_accumulators
  sll_assert                  
  sll_meshes
  sll_memory
  sll_working_precision )

# Ctest
IF(BUILD_TESTING)

#  add_executable(test_picutilities
#                 unit_test.F90 ) 
#  target_link_libraries( test_picutilities
#                         sll_pic_utilities
#                         sll_meshes
#                         sll_assert
#                         sll_memory
#                         sll_working_precision )
#
  add_executable( test_particle_sort_2d unit_test_particle_sort.F90 )
  target_link_libraries( test_particle_sort_2d 
    sll_random_deviate_generators
    sll_pic_utilities
    sll_pic_particles
    sll_meshes )
  ADD_TEST( NAME pic_particle_sort COMMAND test_particle_sort_2d )
  SET_TESTS_PROPERTIES( pic_particle_sort PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(BUILD_TESTING)
