add_library( sll_pic_particles STATIC 
  particle_representation.F90
  particle_group_2d.F90 
  particle_group_4d.F90 )

target_link_libraries( sll_pic_particles
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes )

# Tests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )
