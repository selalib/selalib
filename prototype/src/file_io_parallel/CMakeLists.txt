SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
SET(IO_SRC sll_gnuplot_parallel.F90
           sll_xdmf_serial_blocks.F90)

IF(HDF5_ENABLED AND HDF5_PARALLEL_ENABLED)
   SET(IO_SRC ${IO_SRC} sll_hdf5_io_parallel.F90
                        sll_xdmf_parallel.F90 )
ENDIF()

ADD_LIBRARY(sll_file_io_parallel STATIC ${IO_SRC})

TARGET_LINK_LIBRARIES(sll_file_io_parallel sll_collective
                                           sll_file_io
                                           sll_assertion 
                                           sll_utilities
                                           ${HDF5_LIBRARIES} )


IF(HDF5_ENABLED)
   ADD_EXECUTABLE(test_io_serial_blocks unit_test_serial_blocks.F90)
   TARGET_LINK_LIBRARIES(test_io_serial_blocks sll_file_io_parallel )
ENDIF(HDF5_ENABLED)
