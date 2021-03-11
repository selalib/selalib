add_library(
  sll_meshes STATIC sll_m_meshes_base.F90 sll_m_cartesian_meshes.F90
                    sll_m_hexagonal_meshes.F90 sll_m_triangular_meshes.F90)
target_link_libraries(sll_meshes sll_boundary_condition_descriptors)

target_link_libraries(sll_meshes sll_file_io sll_errors sll_memory sll_assert
                      ${LAPACK_LIBRARIES})

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
