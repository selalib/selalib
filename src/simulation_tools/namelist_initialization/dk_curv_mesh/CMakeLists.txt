add_library(sll_dk_curv_mesh STATIC sll_m_dk_curv_mesh.F90)

target_link_libraries(sll_dk_curv_mesh sll_meshes
                      sll_coordinate_transformations)

# no testing for the moment Ctest IF(BUILD_TESTING)

# ADD_SUBDIRECTORY(testing)

# ENDIF(BUILD_TESTING)
