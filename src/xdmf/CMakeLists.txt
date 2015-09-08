SET( XDMF_SRC sll_xml.F90 sll_xdmf.F90 )
SET( XDMF_DEP sll_working_precision )

IF (MPI_ENABLED AND HDF5_ENABLED AND HDF5_PARALLEL_ENABLED) 
  SET( XDMF_SRC ${XDMF_SRC} sll_xdmf_parallel.F90 )
  SET( XDMF_DEP ${XDMF_DEP} sll_assert sll_collective sll_file_io )
ENDIF()

ADD_LIBRARY( sll_xdmf STATIC ${XDMF_SRC} )
TARGET_LINK_LIBRARIES( sll_xdmf ${XDMF_DEP} )

# Tests
IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( testing )
ENDIF( BUILD_TESTING )

