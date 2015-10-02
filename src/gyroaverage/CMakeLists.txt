SET( GYROAVERAGE_SRC_FILES 
  sll_gyroaverage_2d_base.F90
  sll_gyroaverage_2d_polar.F90
  sll_gyroaverage_utilities.F90
  #sll_gyroaverage_2d_polar_computation.F90
  sll_gyroaverage_2d_polar_pade_solver.F90
  sll_gyroaverage_2d_polar_splines_solver.F90
  sll_gyroaverage_2d_polar_hermite_solver.F90 
  sll_qn_2d_polar_precompute.F90
  sll_qn_2d_base.F90
  sll_qn_2d_polar.F90
  sll_qn_2d_polar_splines_solver.F90)

ADD_LIBRARY( sll_gyroaverage STATIC ${GYROAVERAGE_SRC_FILES} )

TARGET_LINK_LIBRARIES( sll_gyroaverage
  sll_characteristics
  sll_interpolators
  sll_periodic_interpolation
  sll_splines
  sll_file_io
  sll_integration )


# Ctests
IF(BUILD_TESTING)

  SET( INPUT_FILES zeros_bessel.txt )
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

  ADD_EXECUTABLE( test_gyroaverage_2d_polar_hermite
    unit_test_gyroaverage_2d_polar_hermite.F90 )

  ADD_EXECUTABLE( test_gyroaverage_2d_polar_splines
    unit_test_gyroaverage_2d_polar_splines.F90 )

  ADD_EXECUTABLE( test_gyroaverage_2d_polar_pade
    unit_test_gyroaverage_2d_polar_pade.F90 )

  ADD_EXECUTABLE(  test_qn_2d_polar_splines
    unit_test_qn_2d_polar_splines.F90)  
 
  ADD_EXECUTABLE(  test_qn_2d_polar_precompute
    unit_test_qn_2d_polar_precompute.F90)

  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_hermite sll_gyroaverage )
  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_splines sll_gyroaverage )
  TARGET_LINK_LIBRARIES( test_gyroaverage_2d_polar_pade    sll_gyroaverage )
  TARGET_LINK_LIBRARIES( test_qn_2d_polar_splines          sll_gyroaverage )  
  TARGET_LINK_LIBRARIES( test_qn_2d_polar_precompute       sll_gyroaverage )

  ADD_TEST( NAME gyroaverage_polar_hermite COMMAND test_gyroaverage_2d_polar_hermite )
  ADD_TEST( NAME gyroaverage_polar_splines COMMAND test_gyroaverage_2d_polar_splines )
  ADD_TEST( NAME gyroaverage_polar_pade    COMMAND test_gyroaverage_2d_polar_pade    )
  ADD_TEST( NAME qn_polar_splines          COMMAND test_qn_2d_polar_splines          )
  ADD_TEST( NAME qn_polar_precompute       COMMAND test_qn_2d_polar_precompute       )

  SET_TESTS_PROPERTIES( gyroaverage_polar_hermite PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( gyroaverage_polar_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( gyroaverage_polar_pade    PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( qn_polar_splines          PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
  SET_TESTS_PROPERTIES( qn_polar_precompute       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(BUILD_TESTING)



### OLD CMAKELISTS.txt before prototype/src move

#   SET( GYROAVERAGE_SRC_FILES 
# 	sll_gyroaverage_2d_base.F90
# 	sll_gyroaverage_2d_polar.F90
# 	sll_gyroaverage_utilities.F90
# 	#sll_gyroaverage_2d_polar_computation.F90
#         sll_gyroaverage_2d_polar_pade_solver.F90
#         sll_gyroaverage_2d_polar_splines_solver.F90
#         sll_gyroaverage_2d_polar_hermite_solver.F90
#         sll_qn_2d_polar_precompute.F90
#         sll_qn_2d_base.F90
#         sll_qn_2d_polar.F90
#         sll_qn_2d_polar_splines_solver.F90)
# 
#   ADD_LIBRARY(
#     sll_gyroaverage
#     STATIC
#     ${GYROAVERAGE_SRC_FILES})
# 
#   TARGET_LINK_LIBRARIES( 
#   sll_gyroaverage
#   sll_characteristics
#   sll_interpolators
#   sll_periodic_interp
#   sll_splines
#   sll_file_io
#   sll_integration)
# 
#   ADD_EXECUTABLE(
#     test_gyroaverage_2d_polar_hermite
#     unit_test_gyroaverage_2d_polar_hermite.F90)
# 
#   ADD_EXECUTABLE(
#     test_gyroaverage_2d_polar_splines
#     unit_test_gyroaverage_2d_polar_splines.F90)
# 
#   ADD_EXECUTABLE(
#     test_gyroaverage_2d_polar_pade
#     unit_test_gyroaverage_2d_polar_pade.F90)
#     
#   ADD_EXECUTABLE(
#     test_qn_2d_polar_splines
#     unit_test_qn_2d_polar_splines.F90)  
# 
#   ADD_EXECUTABLE(
#     test_qn_2d_polar_precompute
#     unit_test_qn_2d_polar_precompute.F90)  
# 
# 
# 
#   TARGET_LINK_LIBRARIES(
#     test_gyroaverage_2d_polar_hermite
#     sll_gyroaverage)
# 
#   TARGET_LINK_LIBRARIES(
#     test_gyroaverage_2d_polar_splines
#     sll_gyroaverage)
# 
#   TARGET_LINK_LIBRARIES(
#     test_gyroaverage_2d_polar_pade
#     sll_gyroaverage)
#     
#   TARGET_LINK_LIBRARIES(
#     test_qn_2d_polar_splines
#     sll_gyroaverage)  
# 
#   TARGET_LINK_LIBRARIES(
#     test_qn_2d_polar_precompute
#     sll_gyroaverage)  
# 


