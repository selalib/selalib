IF(HDF5_IS_PARALLEL)
   MESSAGE(STATUS "Switch to MPI fortran Wrapper to build sll_file_io library")
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ELSE()
   MESSAGE(STATUS "Serial version of HDF5 is installed")
ENDIF(HDF5_IS_PARALLEL)

SET(IO_SRC sll_ascii_io.F90 
           sll_binary_io.F90 
           sll_gnuplot.F90 
           sll_hdf5_io_serial.F90
           sll_xml_io.F90
           sll_plotmtv.F90
           sll_tri_mesh_xmf.F90
           sll_xdmf.F90 )

ADD_LIBRARY(sll_file_io STATIC ${IO_SRC})

TARGET_LINK_LIBRARIES(sll_file_io sll_utilities ${HDF5_LIBRARIES})

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
