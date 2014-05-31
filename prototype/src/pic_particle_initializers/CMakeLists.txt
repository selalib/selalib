#HACK: THIS SHOULD BE MADE BETTER AND AT A LOWER LEVEL. Problem: when hdf5 is
# present and is parallel, we need to use the MPI compiler wrapper and this
# should propagate to all its clients, even if these were envisioned as
# sequential modules. Or, maybe this should be handled at the highest possible
# level: once the parallel hdf5 is found, all the library should be compiled 
# with the mpi wrapper compiler.
IF(HDF5_IS_PARALLEL)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)

ADD_LIBRARY( sll_pic_initializers STATIC 
  sll_particle_init4D.F90 )

TARGET_LINK_LIBRARIES( sll_pic_initializers
  sll_pic_particles
  sll_visu_pic
  sll_logical_meshes
  sll_working_precision
  sll_assertion
  sll_memory )

ADD_EXECUTABLE(test_pic_initializers unit_test.F90  )

TARGET_LINK_LIBRARIES(test_pic_initializers
  sll_pic_initializers
#  sll_visu_pic
 sll_collective
 # sll_file_io )
)
