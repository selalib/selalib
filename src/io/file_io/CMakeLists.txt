IF(HDF5_IS_PARALLEL)
   MESSAGE(STATUS "Switch to MPI fortran Wrapper to build sll_file_io library")
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ELSE()
   MESSAGE(STATUS "Serial version of HDF5 is installed")
ENDIF(HDF5_IS_PARALLEL)

SET(IO_SRC sll_m_ascii_io.F90 
           sll_m_binary_io.F90 
           sll_m_gnuplot.F90 
           sll_m_hdf5_io_serial.F90
           sll_m_xml_io.F90
           sll_m_plotmtv.F90
           sll_m_tri_mesh_xmf.F90
           sll_m_xdmf.F90 )

ADD_LIBRARY(sll_file_io STATIC ${IO_SRC})

TARGET_INCLUDE_DIRECTORIES(sll_file_io BEFORE PUBLIC ${HDF5_INCLUDE_DIR_FORTRAN})
TARGET_LINK_LIBRARIES(sll_file_io sll_utilities ${HDF5_LIBRARIES})

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
