IF(HDF5_IS_PARALLEL AND HDF5_ENABLED)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL AND HDF5_ENABLED)

ADD_EXECUTABLE( test_pic_visu test_pic_visu.F90 biot_savart.F90 )
TARGET_LINK_LIBRARIES( test_pic_visu sll_pic_visu )
SET( PROCS 1 )
SET( ARGS " ")
ADD_MPI_TEST( pic_visu test_pic_visu ${PROCS} ${ARGS} )
SET( pic_visu PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_pic_viewer test_pic_viewer.F90 biot_savart.F90 )
TARGET_LINK_LIBRARIES( test_pic_viewer sll_pic_visu 
                                       sll_poisson_solvers
                                       sll_file_io
                                       sll_coordinate_transformations )
SET( PROCS 1 )
SET( ARGS " ")
ADD_MPI_TEST( pic_viewer test_pic_viewer ${PROCS} ${ARGS} )
SET( pic_visu PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
