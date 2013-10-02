MESSAGE(STATUS "HDF5_IS_PARALLEL:${HDF5_IS_PARALLEL}")
IF(HDF5_IS_PARALLEL)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)

SET(IO_SRC sll_ascii_io.F90 
           sll_binary_io.F90 
           sll_gnuplot.F90 
           sll_hdf5_io.F90
           sll_xml_io.F90
           sll_plotmtv.F90
           sll_xdmf.F90 )

ADD_LIBRARY(sll_file_io STATIC ${IO_SRC})

ADD_EXECUTABLE(test_io unit_test.F90)

TARGET_LINK_LIBRARIES(sll_file_io sll_utilities ${HDF5_LIBRARIES})

TARGET_LINK_LIBRARIES(test_io sll_file_io sll_constants)

ADD_TEST(NAME file_io COMMAND test_io)
SET_TESTS_PROPERTIES(file_io PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

IF(HDF5_PARALLEL_ENABLED)

   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   ADD_DEFINITIONS(-DHDF5_PARALLEL)
   SET(IO_SRC ${IO_SRC} sll_hdf5_io_parallel.F90
                        sll_gnuplot_parallel.F90
                        sll_xdmf_parallel.F90 )

   ADD_LIBRARY(sll_file_io_parallel STATIC ${IO_SRC})

   TARGET_LINK_LIBRARIES(sll_file_io_parallel sll_utilities 
                                              sll_collective
                                              ${HDF5_LIBRARIES} )
   IF(NOT STDF95) 

      ADD_EXECUTABLE(test_io_parallel unit_test_parallel.F90)

      TARGET_LINK_LIBRARIES(test_io_parallel sll_file_io_parallel
                                             sll_remap)
   ENDIF()

ENDIF()
