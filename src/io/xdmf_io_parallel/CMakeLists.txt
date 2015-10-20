SET( LIB sll_xdmf_io_parallel )
SET( SRC sll_xdmf_io_parallel.F90 sll_m_xdmf_light_parallel.F90 )
SET( DEP sll_xdmf_io sll_collective )

ADD_LIBRARY( ${LIB} STATIC ${SRC} )
TARGET_LINK_LIBRARIES( ${LIB} ${DEP} )

# Tests
IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( testing )
ENDIF( BUILD_TESTING )

