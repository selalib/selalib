set(SRC_FILES sll_m_nml_mesh_1d_landau_cart.F90 sll_m_nml_mesh_1d_unif_cart.F90
              sll_m_nml_mesh_1d_two_grid_cart.F90 sll_m_nml_mesh_1d_cart.F90)

add_library(sll_nml_mesh_cart STATIC ${SRC_FILES})

target_link_libraries(sll_nml_mesh_cart sll_meshes)

# Ctest
if(BUILD_TESTING)

  add_subdirectory(testing)

endif(BUILD_TESTING)
