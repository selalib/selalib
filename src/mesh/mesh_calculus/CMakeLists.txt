add_library(sll_mesh_calculus STATIC sll_m_mesh_calculus_2d.F90)

target_link_libraries(sll_mesh_calculus sll_coordinate_transformations
                      sll_meshes sll_quadrature sll_assert)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
