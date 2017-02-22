SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
SET(IO_SRC sll_m_gnuplot_parallel.F90 sll_m_xdmf_serial_blocks.F90)

IF(HDF5_ENABLED AND HDF5_PARALLEL_ENABLED)
   SET(IO_SRC ${IO_SRC} sll_m_hdf5_io_parallel.F90
                        sll_m_xdmf_parallel.F90 )
ENDIF()

ADD_LIBRARY(sll_file_io_parallel STATIC ${IO_SRC})

TARGET_LINK_LIBRARIES(sll_file_io_parallel sll_collective
                                           sll_file_io
                                           sll_assert 
                                           sll_utilities
                                           ${HDF5_LIBRARIES} )

IF(HDF5_ENABLED AND BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(HDF5_ENABLED AND BUILD_TESTING)
