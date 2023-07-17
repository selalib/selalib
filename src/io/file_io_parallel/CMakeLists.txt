set(IO_SRC sll_m_gnuplot_parallel.F90 sll_m_xdmf_serial_blocks.F90)

if(HDF5_ENABLED AND HDF5_PARALLEL_ENABLED)
  set(IO_SRC ${IO_SRC} sll_m_hdf5_io_parallel.F90 sll_m_xdmf_parallel.F90)
endif()

add_library(sll_file_io_parallel STATIC ${IO_SRC})

target_include_directories(sll_file_io_parallel BEFORE
                           PUBLIC ${HDF5_INCLUDE_DIR_FORTRAN})

target_link_libraries(sll_file_io_parallel sll_collective sll_file_io
                      sll_assert sll_utilities ${HDF5_LIBRARIES} ${HDF5_Fortran_LIBRARIES})

if(HDF5_ENABLED AND BUILD_TESTING)
  add_subdirectory(testing)
endif(HDF5_ENABLED AND BUILD_TESTING)
