#HACK: THIS SHOULD BE MADE BETTER AND AT A LOWER LEVEL. Problem: when hdf5 is
# present and is parallel, we need to use the MPI compiler wrapper and this
# should propagate to all its clients, even if these were envisioned as
# sequential modules. visu_pic is an example.

IF(HDF5_IS_PARALLEL AND HDF5_ENABLED)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL AND HDF5_ENABLED)

ADD_LIBRARY(sll_visu_pic STATIC sll_visu_pic.F90)

TARGET_LINK_LIBRARIES(sll_visu_pic sll_utilities
                                   sll_file_io
                                   sll_assertion
                                   sll_working_precision )

ADD_EXECUTABLE(test_visu_pic unit_test.F90 biot_savart.F90)

TARGET_LINK_LIBRARIES(test_visu_pic sll_visu_pic)
