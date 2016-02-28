SET(MEMORY_SRC sll_m_memory.F90)
ADD_LIBRARY(sll_memory STATIC ${MEMORY_SRC})
TARGET_LINK_LIBRARIES(sll_memory sll_working_precision sll_constants)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
