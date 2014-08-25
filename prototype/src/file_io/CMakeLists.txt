IF(HDF5_ENABLED AND HDF5_IS_PARALLEL)
   MESSAGE(STATUS "HDF5_IS_PARALLEL:${HDF5_IS_PARALLEL}")
   MESSAGE(STATUS "Switch to MPI fortran Wrapper to build sll_file_io library")
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL)

SET(IO_SRC sll_ascii_io.F90 
           sll_binary_io.F90 
           sll_gnuplot.F90 
           sll_hdf5_io_serial.F90
           sll_xml_io.F90
           sll_plotmtv.F90
           sll_xdmf.F90 )

ADD_LIBRARY(sll_file_io STATIC ${IO_SRC})

ADD_EXECUTABLE(test_io unit_test.F90)

TARGET_LINK_LIBRARIES(sll_file_io sll_utilities ${HDF5_LIBRARIES})

TARGET_LINK_LIBRARIES(test_io sll_file_io sll_constants)

ADD_TEST(NAME file_io COMMAND test_io)

SET_TESTS_PROPERTIES(file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_io_tri_mesh unit_test_tri_mesh.F90 
                                sll_tri_mesh_xmf.F90
                                sll_generate_tri_mesh.F90)
