SET( GYROAVERAGE_SRC_FILES 
  sll_m_gyroaverage_2d_base.F90
  sll_m_gyroaverage_2d_polar.F90
  sll_m_gyroaverage_utilities.F90
  #sll_m_gyroaverage_2d_polar_computation.F90
  sll_m_gyroaverage_2d_polar_pade_solver.F90
  sll_m_gyroaverage_2d_polar_splines_solver.F90
  sll_m_gyroaverage_2d_polar_hermite_solver.F90 
  sll_m_qn_2d_polar_precompute.F90
  sll_m_qn_2d_base.F90
  sll_m_qn_2d_polar.F90
  sll_m_qn_2d_polar_splines_solver.F90)

ADD_LIBRARY( sll_gyroaverage STATIC ${GYROAVERAGE_SRC_FILES} )

TARGET_LINK_LIBRARIES( sll_gyroaverage
  sll_characteristics
  sll_interpolators
  sll_periodic_interpolation
  sll_splines
  sll_file_io
  sll_quadrature 
  ${BLAS_LIBRARIES})


# Ctests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)



### OLD CMAKELISTS.txt before prototype/src move

#   SET( GYROAVERAGE_SRC_FILES 
# 	sll_m_gyroaverage_2d_base.F90
# 	sll_m_gyroaverage_2d_polar.F90
# 	sll_m_gyroaverage_utilities.F90
# 	#sll_m_gyroaverage_2d_polar_computation.F90
#         sll_m_gyroaverage_2d_polar_pade_solver.F90
#         sll_m_gyroaverage_2d_polar_splines_solver.F90
#         sll_m_gyroaverage_2d_polar_hermite_solver.F90
#         sll_m_qn_2d_polar_precompute.F90
#         sll_m_qn_2d_base.F90
#         sll_m_qn_2d_polar.F90
#         sll_m_qn_2d_polar_splines_solver.F90)
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


