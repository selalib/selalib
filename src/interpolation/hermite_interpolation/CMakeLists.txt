add_library(
  sll_hermite_interpolation STATIC
  sll_m_hermite_interpolation_1d.F90 sll_m_hermite_interpolation_2d.F90
  sll_m_interpolation_hex_hermite.F90
  sll_m_hermite_aligned_interpolation_2d.F90)

target_link_libraries(sll_hermite_interpolation sll_meshes sll_utilities)

# Ctests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
