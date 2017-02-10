ADD_EXECUTABLE( test_advection_1d_spectral test_advection_1d_spectral.F90)
TARGET_LINK_LIBRARIES( test_advection_1d_spectral sll_advection sll_file_io sll_fft)

ADD_EXECUTABLE( test_advection_1d_ampere test_advection_1d_ampere.F90)
TARGET_LINK_LIBRARIES( test_advection_1d_ampere sll_advection sll_file_io sll_fft)

ADD_EXECUTABLE( test_advection_1d_periodic test_advection_1d_periodic.F90)
TARGET_LINK_LIBRARIES( test_advection_1d_periodic sll_advection)
ADD_TEST( NAME advection_1d_periodic COMMAND test_advection_1d_periodic )
SET_TESTS_PROPERTIES( advection_1d_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_advection_1d_non_uniform_cubic_splines test_advection_1d_non_uniform_cubic_splines.F90)
TARGET_LINK_LIBRARIES( test_advection_1d_non_uniform_cubic_splines sll_advection)
ADD_TEST( NAME advection_1d_non_uniform_cubic_splines COMMAND test_advection_1d_non_uniform_cubic_splines )
SET_TESTS_PROPERTIES( advection_1d_non_uniform_cubic_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_advection_1d_bsl test_advection_1d_bsl.F90)
TARGET_LINK_LIBRARIES( test_advection_1d_bsl sll_advection)
#ADD_TEST( NAME advection_1d_bsl COMMAND test_advection_1d_bsl      )
#SET_TESTS_PROPERTIES( advection_1d_bsl  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_advection_1d_csl_periodic test_advection_1d_csl_periodic.F90)
TARGET_LINK_LIBRARIES( test_advection_1d_csl_periodic sll_advection)

ADD_EXECUTABLE( test_advection_1d_csl test_advection_1d_csl.F90)
TARGET_LINK_LIBRARIES( test_advection_1d_csl sll_advection)
#ADD_TEST( NAME advection_1d_csl COMMAND test_advection_1d_csl      )
#SET_TESTS_PROPERTIES( advection_1d_csl PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" ) # ISSUE

ADD_EXECUTABLE( test_advection_1d_psm test_advection_1d_psm.F90)
TARGET_LINK_LIBRARIES( test_advection_1d_psm sll_advection)
#ADD_TEST( NAME advection_1d_psm COMMAND test_advection_1d_psm      )
#SET_TESTS_PROPERTIES( advection_1d_psm PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" ) # ISSUE

ADD_EXECUTABLE( test_advection_2d_bsl test_advection_2d_bsl.F90)
TARGET_LINK_LIBRARIES( test_advection_2d_bsl sll_advection)
ADD_TEST( NAME advection_2d_bsl COMMAND test_advection_2d_bsl      )
SET_TESTS_PROPERTIES( advection_2d_bsl  PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_advection_2d_integer_oblic test_advection_2d_integer_oblic.F90)
TARGET_LINK_LIBRARIES( test_advection_2d_integer_oblic sll_advection)

ADD_EXECUTABLE( test_advection_2d_oblic test_advection_2d_oblic.F90)
TARGET_LINK_LIBRARIES( test_advection_2d_oblic sll_advection)

ADD_EXECUTABLE( test_advection_2d_tensor_product test_advection_2d_tensor_product.F90)
TARGET_LINK_LIBRARIES( test_advection_2d_tensor_product sll_advection)
ADD_TEST( NAME advection_2d_tensor_product COMMAND test_advection_2d_tensor_product )
SET_TESTS_PROPERTIES( advection_2d_tensor_product PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_advection_2d_csl test_advection_2d_csl.F90)
TARGET_LINK_LIBRARIES( test_advection_2d_csl sll_advection)
#ADD_TEST( NAME advection_2d_csl COMMAND test_advection_2d_csl      )
#SET_TESTS_PROPERTIES( advection_2d_csl PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" ) # ISSUE

ADD_EXECUTABLE(test_advection_2d_tri_mesh test_advection_2d_tri_mesh.F90)
TARGET_LINK_LIBRARIES(test_advection_2d_tri_mesh sll_advection sll_mesh_calculus)

IF(MPI_ENABLED AND HDF5_ENABLED)

  ADD_EXECUTABLE( test_aligned_translation_2d aligned_translation_2d.F90)

  TARGET_LINK_LIBRARIES( test_aligned_translation_2d sll_file_io sll_advection
       sll_utilities sll_fcisl)

ENDIF(MPI_ENABLED AND HDF5_ENABLED)

ADD_EXECUTABLE( test_rotation_2d rotation_2d.F90)
TARGET_LINK_LIBRARIES( test_rotation_2d sll_advection)

ADD_EXECUTABLE( test_advection_2d_nufft test_advection_2d_nufft.F90)
TARGET_LINK_LIBRARIES( test_advection_2d_nufft sll_advection)
