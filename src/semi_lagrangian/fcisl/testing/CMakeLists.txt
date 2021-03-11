add_executable(test_fcisl test_fcisl.F90)
target_link_libraries(test_fcisl sll_fcisl sll_meshes sll_advection)

add_executable(test_derivative_2d_oblic test_derivative_2d_oblic.F90)
target_link_libraries(test_derivative_2d_oblic sll_fcisl sll_meshes
                      sll_advection)

add_executable(test_fcisl_toroidal test_fcisl_toroidal.F90)
target_link_libraries(test_fcisl_toroidal sll_fcisl sll_meshes sll_advection)

add_executable(test_aligned_derivative_2d aligned_derivative_2d.F90)
target_link_libraries(test_aligned_derivative_2d sll_fcisl sll_meshes
                      sll_advection)
