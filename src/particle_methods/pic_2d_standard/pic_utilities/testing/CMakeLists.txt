# CTest
add_executable( test_particle_sort_2d unit_test_particle_sort.F90 )
target_link_libraries( test_particle_sort_2d 
  sll_random_deviate_generators
  sll_pic_utilities
  sll_pic_particles
  sll_meshes )
add_test( NAME pic_particle_sort COMMAND test_particle_sort_2d )
set_tests_properties( pic_particle_sort PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
