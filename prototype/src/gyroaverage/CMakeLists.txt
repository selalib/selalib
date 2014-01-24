IF( NOT STDF95)

  SET( GYROAVERAGE_SRC_FILES 
	sll_gyroaverage_2d_polar_computation.F90
	sll_gyroaverage_2d_polar.F90
	sll_gyroaverage_2d_base.F90)


  ADD_LIBRARY(
    sll_gyroaverage
    STATIC
    ${GYROAVERAGE_SRC_FILES})

  TARGET_LINK_LIBRARIES( 
  sll_gyroaverage
  sll_periodic_interp
  sll_splines
  sll_characteristics
  sll_interpolators
)

  ADD_EXECUTABLE(
    test_gyroaverage_2d_polar_computation
    unit_test_gyroaverage_2d_polar_computation.F90)

  TARGET_LINK_LIBRARIES(
    test_gyroaverage_2d_polar_computation
    sll_gyroaverage)




ENDIF(NOT STDF95)
