SET( SRC_FILES    
  sll_m_nml_mesh_1d_landau_cart.F90
  sll_m_nml_mesh_1d_unif_cart.F90
  sll_m_nml_mesh_1d_two_grid_cart.F90
  sll_m_nml_mesh_1d_cart.F90 )

ADD_LIBRARY( sll_nml_mesh_cart STATIC
             ${SRC_FILES})

TARGET_LINK_LIBRARIES( sll_nml_mesh_cart sll_meshes )

# Ctest
IF(BUILD_TESTING)

  ADD_SUBDIRECTORY(testing)


ENDIF(BUILD_TESTING)
