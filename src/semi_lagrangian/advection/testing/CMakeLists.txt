add_executable(test_advection_1d_spectral test_advection_1d_spectral.F90)
target_link_libraries(test_advection_1d_spectral sll_advection sll_file_io
                      sll_fft)

add_executable(test_advection_1d_ampere test_advection_1d_ampere.F90)
target_link_libraries(test_advection_1d_ampere sll_advection sll_file_io
                      sll_fft)

add_executable(test_advection_1d_periodic test_advection_1d_periodic.F90)
target_link_libraries(test_advection_1d_periodic sll_advection)
add_test(NAME advection_1d_periodic COMMAND test_advection_1d_periodic)
set_tests_properties(advection_1d_periodic PROPERTIES PASS_REGULAR_EXPRESSION
                                                      "PASSED")

add_executable(test_advection_1d_non_uniform_cubic_splines
               test_advection_1d_non_uniform_cubic_splines.F90)
target_link_libraries(test_advection_1d_non_uniform_cubic_splines sll_advection)
add_test(NAME advection_1d_non_uniform_cubic_splines
         COMMAND test_advection_1d_non_uniform_cubic_splines)
set_tests_properties(advection_1d_non_uniform_cubic_splines
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_advection_1d_bsl test_advection_1d_bsl.F90)
target_link_libraries(test_advection_1d_bsl sll_advection)
# ADD_TEST( NAME advection_1d_bsl COMMAND test_advection_1d_bsl      )
# SET_TESTS_PROPERTIES( advection_1d_bsl  PROPERTIES PASS_REGULAR_EXPRESSION
# "PASSED" )

add_executable(test_advection_1d_csl_periodic
               test_advection_1d_csl_periodic.F90)
target_link_libraries(test_advection_1d_csl_periodic sll_advection)

add_executable(test_advection_1d_csl test_advection_1d_csl.F90)
target_link_libraries(test_advection_1d_csl sll_advection)
# ADD_TEST( NAME advection_1d_csl COMMAND test_advection_1d_csl      )
# SET_TESTS_PROPERTIES( advection_1d_csl PROPERTIES PASS_REGULAR_EXPRESSION
# "PASSED" ) # ISSUE

add_executable(test_advection_1d_psm test_advection_1d_psm.F90)
target_link_libraries(test_advection_1d_psm sll_advection)
# ADD_TEST( NAME advection_1d_psm COMMAND test_advection_1d_psm      )
# SET_TESTS_PROPERTIES( advection_1d_psm PROPERTIES PASS_REGULAR_EXPRESSION
# "PASSED" ) # ISSUE

add_executable(test_advection_2d_bsl test_advection_2d_bsl.F90)
target_link_libraries(test_advection_2d_bsl sll_advection)
add_test(NAME advection_2d_bsl COMMAND test_advection_2d_bsl)
set_tests_properties(advection_2d_bsl PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED")

add_executable(test_advection_2d_integer_oblic
               test_advection_2d_integer_oblic.F90)
target_link_libraries(test_advection_2d_integer_oblic sll_advection)

add_executable(test_advection_2d_oblic test_advection_2d_oblic.F90)
target_link_libraries(test_advection_2d_oblic sll_advection)

add_executable(test_advection_2d_tensor_product
               test_advection_2d_tensor_product.F90)
target_link_libraries(test_advection_2d_tensor_product sll_advection)
add_test(NAME advection_2d_tensor_product
         COMMAND test_advection_2d_tensor_product)
set_tests_properties(advection_2d_tensor_product
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_advection_2d_csl test_advection_2d_csl.F90)
target_link_libraries(test_advection_2d_csl sll_advection)
# ADD_TEST( NAME advection_2d_csl COMMAND test_advection_2d_csl      )
# SET_TESTS_PROPERTIES( advection_2d_csl PROPERTIES PASS_REGULAR_EXPRESSION
# "PASSED" ) # ISSUE

add_executable(test_advection_2d_tri_mesh test_advection_2d_tri_mesh.F90)
target_link_libraries(test_advection_2d_tri_mesh sll_advection
                      sll_mesh_calculus)

if(MPI_ENABLED
   AND HDF5_ENABLED
   AND HDF5_IS_PARALLEL)

  add_executable(test_aligned_translation_2d aligned_translation_2d.F90)

  target_link_libraries(test_aligned_translation_2d sll_file_io sll_advection
                        sll_utilities sll_fcisl)

endif(
  MPI_ENABLED
  AND HDF5_ENABLED
  AND HDF5_IS_PARALLEL)

add_executable(test_rotation_2d rotation_2d.F90)
target_link_libraries(test_rotation_2d sll_advection)

add_executable(test_advection_2d_nufft test_advection_2d_nufft.F90)
target_link_libraries(test_advection_2d_nufft sll_advection)
