add_library(
  sll_lagrange_interpolation STATIC
  sll_m_lagrange_interpolation.F90 sll_m_lagrange_interpolation_1d.F90
  sll_m_lagrange_interpolation_1d_fast.F90)

target_link_libraries(sll_lagrange_interpolation sll_utilities
                      sll_boundary_condition_descriptors)

# Ctest
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
