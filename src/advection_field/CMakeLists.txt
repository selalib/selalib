#ADD_LIBRARY(sll_advection_field STATIC sll_advection_field.F90)

#TARGET_LINK_LIBRARIES(sll_advection_field sll_distribution_function)

# Ctest
#IF(BUILD_TESTING)
#  ADD_EXECUTABLE( test_advection_field unit_test.F90 )
#  TARGET_LINK_LIBRARIES( test_advection_field sll_advection_field )
#  ADD_TEST( NAME advection_field COMMAND test_advection_field )
#ENDIF(BUILD_TESTING)
