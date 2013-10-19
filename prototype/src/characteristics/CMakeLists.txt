IF( NOT STDF95)

  SET( CHARACTERISTICS_SRC_FILES 
	sll_characteristics_2d_explicit_euler.F90
	sll_characteristics_2d_base.F90)


  ADD_LIBRARY(
    sll_characteristics
    STATIC
    ${CHARACTERISTICS_SRC_FILES})

  TARGET_LINK_LIBRARIES( 
    sll_characteristics
    sll_utilities)

  ADD_EXECUTABLE(
    test_characteristics
    unit_test.F90)

  TARGET_LINK_LIBRARIES(
    test_characteristics
    sll_characteristics)



ENDIF(NOT STDF95)
