add_library(sll_decomposition STATIC sll_m_decomposition.F90 sll_m_decomposition_advanced.F90)
IF(USE_FMEMPOOL)
  set( FMEMPOOL "fmempool" )
  add_definitions( -DUSE_FMEMPOOL )
ENDIF(USE_FMEMPOOL)

target_link_libraries(sll_decomposition sll_collective ${FMEMPOOL} sll_compression )

# Ctests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
