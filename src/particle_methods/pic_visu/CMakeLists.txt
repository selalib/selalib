#HACK: THIS SHOULD BE MADE BETTER AND AT A LOWER LEVEL. Problem: when hdf5 is
# present and is parallel, we need to use the MPI compiler wrapper and this
# should propagate to all its clients, even if these were envisioned as
# sequential modules. pic_visu is an example.

IF(HDF5_IS_PARALLEL AND HDF5_ENABLED)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL AND HDF5_ENABLED)

ADD_LIBRARY( sll_pic_visu STATIC sll_m_pic_visu.F90 )
TARGET_LINK_LIBRARIES( sll_pic_visu
  sll_utilities
  sll_file_io
  sll_assert)

    
IF(BUILD_TESTING AND MPI_ENABLED)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING AND MPI_ENABLED)
