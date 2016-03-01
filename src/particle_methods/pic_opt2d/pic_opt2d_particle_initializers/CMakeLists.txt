#HACK: THIS SHOULD BE MADE BETTER AND AT A LOWER LEVEL. Problem: when hdf5 is
# present and is parallel, we need to use the MPI compiler wrapper and this
# should propagate to all its clients, even if these were envisioned as
# sequential modules. Or, maybe this should be handled at the highest possible
# level: once the parallel hdf5 is found, all the library should be compiled 
# with the mpi wrapper compiler.
if( HDF5_ENABLED AND HDF5_IS_PARALLEL )
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
endif( HDF5_ENABLED AND HDF5_IS_PARALLEL )

add_library( sll_pic_opt2d_initializers STATIC 
  sll_m_particle_initializers_4d.F90 
  sll_m_particle_initializers_2d.F90 )

target_link_libraries( sll_pic_opt2d_initializers
  sll_pic_opt2d_particles
  sll_meshes
  sll_assert
  sll_memory
  sll_random_deviate_generators )

# Tests
if( MPI_ENABLED AND BUILD_TESTING )
  add_subdirectory( testing )
endif( MPI_ENABLED AND BUILD_TESTING )
