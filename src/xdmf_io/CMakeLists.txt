SET( XDMF_SRC sll_m_xdmf.F90 sll_m_xml.F90 sll_m_xdmf_sequential.F90 sll_m_io_utilities.F90 )
SET( XDMF_DEP sll_working_precision )

# Object-oriented wrapper around Pierre's "sll_hdf5_io_serial"
SET( XDMF_SRC ${XDMF_SRC} sll_m_hdf5_serial.F90 )
SET( XDMF_DEP ${XDMF_DEP} sll_file_io )

IF (MPI_ENABLED AND HDF5_ENABLED AND HDF5_PARALLEL_ENABLED)
  SET( XDMF_SRC ${XDMF_SRC} sll_m_xdmf_parallel.F90 )
  SET( XDMF_DEP ${XDMF_DEP} sll_collective )
ENDIF()

ADD_LIBRARY( sll_xdmf STATIC ${XDMF_SRC} )
TARGET_LINK_LIBRARIES( sll_xdmf ${XDMF_DEP} )

# Tests
IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( testing )
ENDIF( BUILD_TESTING )

