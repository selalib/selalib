# CTest
add_executable( test_particle_sort test_particle_sort.F90 )
target_link_libraries( test_particle_sort
  sll_random_deviate_generators
  sll_pic_opt2d_utilities
  sll_pic_opt2d_particles
  sll_meshes )
add_test( NAME pic_particle_sort COMMAND test_particle_sort )
set_tests_properties( pic_particle_sort PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
