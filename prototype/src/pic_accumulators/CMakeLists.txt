#HACK: THIS SHOULD BE MADE BETTER AND AT A LOWER LEVEL. Problem: when hdf5 is
# present and is parallel, we need to use the MPI compiler wrapper and this
# should propagate to all its clients, even if these were envisioned as
# sequential modules. Or, maybe this should be handled at the highest possible
# level: once the parallel hdf5 is found, all the library should be compiled 
# with the mpi wrapper compiler.
IF(HDF5_IS_PARALLEL)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)



ADD_LIBRARY( sll_pic_accumulators STATIC 
  sll_accumulators.F90 )

TARGET_LINK_LIBRARIES(sll_pic_accumulators
  sll_logical_meshes
  sll_memory
  sll_assertion
  sll_working_precision )


ADD_EXECUTABLE(test_pic_accumulator unit_test.F90  )

TARGET_LINK_LIBRARIES(test_pic_accumulator
  sll_pic_utilities
  sll_pic_initializers
  sll_pic_accumulators
  sll_pic_particles
  sll_logical_meshes
  sll_assertion )
#  sll_pic_particles
#  sll_pic_initializers )

