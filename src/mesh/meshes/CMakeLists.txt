ADD_LIBRARY( sll_meshes STATIC sll_m_meshes_base.F90
                               sll_m_cartesian_meshes.F90
                               sll_m_hexagonal_meshes.F90
                               sll_m_triangular_meshes.F90 )
TARGET_LINK_LIBRARIES(sll_meshes sll_boundary_condition_descriptors)

TARGET_LINK_LIBRARIES( sll_meshes sll_file_io
                                  sll_memory
                                  sll_assert )

# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
