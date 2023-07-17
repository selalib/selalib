if(NOT HDF5_IS_PARALLEL)
  message(STATUS "Serial version of HDF5 is installed")
endif(NOT HDF5_IS_PARALLEL)

set(IO_SRC
    sll_m_ascii_io.F90
    sll_m_binary_io.F90
    sll_m_gnuplot.F90
    sll_m_hdf5_io_serial.F90
    sll_m_xml_io.F90
    sll_m_plotmtv.F90
    sll_m_tri_mesh_xmf.F90
    sll_m_xdmf.F90)

add_library(sll_file_io STATIC ${IO_SRC})

target_include_directories(sll_file_io BEFORE
                           PUBLIC ${HDF5_Fortran_INCLUDE_DIRS})
target_link_libraries(sll_file_io sll_utilities ${HDF5_LIBRARIES} ${HDF5_Fortran_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
