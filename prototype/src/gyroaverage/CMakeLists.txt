IF( NOT STDF95)

  SET( GYROAVERAGE_SRC_FILES 
	sll_gyroaverage_2d_base.F90
	sll_gyroaverage_2d_polar.F90
	#sll_gyroaverage_2d_polar_computation.F90
        sll_gyroaverage_2d_polar_pade_solver.F90
        sll_gyroaverage_2d_polar_splines_solver.F90
        sll_gyroaverage_2d_polar_hermite_solver.F90)

  ADD_LIBRARY(
    sll_gyroaverage
    STATIC
    ${GYROAVERAGE_SRC_FILES})

  TARGET_LINK_LIBRARIES( 
  sll_gyroaverage
  sll_characteristics
  sll_interpolators
  sll_periodic_interp
  sll_splines)

  ADD_EXECUTABLE(
    test_gyroaverage_2d_polar_hermite
    unit_test_gyroaverage_2d_polar_hermite.F90)

  ADD_EXECUTABLE(
    test_gyroaverage_2d_polar_splines
    unit_test_gyroaverage_2d_polar_splines.F90)

  ADD_EXECUTABLE(
    test_gyroaverage_2d_polar_pade
    unit_test_gyroaverage_2d_polar_pade.F90)

  TARGET_LINK_LIBRARIES(
    test_gyroaverage_2d_polar_hermite
    sll_gyroaverage)

  TARGET_LINK_LIBRARIES(
    test_gyroaverage_2d_polar_splines
    sll_gyroaverage)

  TARGET_LINK_LIBRARIES(
    test_gyroaverage_2d_polar_pade
    sll_gyroaverage)


ENDIF(NOT STDF95)
