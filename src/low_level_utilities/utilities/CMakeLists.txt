SET( SRC_UTILITIES
  sll_m_utilities.F90
  sll_m_fornberg.F90
)

ADD_LIBRARY( sll_utilities STATIC ${SRC_UTILITIES} )
TARGET_LINK_LIBRARIES( sll_utilities sll_errors sll_memory sll_assert )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
