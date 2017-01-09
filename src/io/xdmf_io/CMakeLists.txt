SET( XDMF_SRC sll_xdmf_io.F90 sll_m_xml.F90 sll_m_xdmf_light_serial.F90 sll_m_io_utilities.F90 )

# Object-oriented wrapper around Pierre's "sll_hdf5_io_serial"
SET( XDMF_SRC ${XDMF_SRC} sll_m_hdf5_serial.F90 )
SET( XDMF_DEP sll_file_io )

ADD_LIBRARY( sll_xdmf_io STATIC ${XDMF_SRC} )
TARGET_INCLUDE_DIRECTORIES(sll_xdmf_io BEFORE PRIVATE ${HDF5_INCLUDE_DIR_FORTRAN})
TARGET_LINK_LIBRARIES( sll_xdmf_io ${XDMF_DEP} )

# Tests
IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( testing )
ENDIF( BUILD_TESTING )

