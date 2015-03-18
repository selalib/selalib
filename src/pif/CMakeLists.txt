ADD_LIBRARY( sll_pif_fieldsolver STATIC sll_pif_fieldsolver.F90 )
TARGET_LINK_LIBRARIES( sll_pif_fieldsolver
  sll_memory
  sll_assert
  sll_working_precision
  sll_collective
  sll_constants
  sll_utilities)

ADD_LIBRARY( sll_wedge_product_generaldim STATIC sll_wedge_product_generaldim.F90 )
TARGET_LINK_LIBRARIES( sll_wedge_product_generaldim
  sll_memory
  sll_assert
  sll_working_precision)
  
  
ADD_EXECUTABLE(pif_fieldsolver_test  pif_fieldsolver_test.F90)
TARGET_LINK_LIBRARIES(pif_fieldsolver_test sll_pif_fieldsolver
					      sll_timer)

# ADD_EXECUTABLE(sll_pif_general  sll_pif_general.F90)
# TARGET_LINK_LIBRARIES(sll_pif_general sll_pif_fieldsolver
# 					      sll_timer 
# 					      sll_memory
# 					      sll_utilities
# 					      sll_constants
# 					      sll_sobol
# 					      sll_prob
# 					      sll_visu_pic
# 					      sll_collective
# 					      sll_moment_matching
# 					      sll_pic_utilities)

					      
# Ctest
# IF(BUILD_TESTING)
#   ADD_EXECUTABLE( test_pic_accumulator unit_test.F90 )
#   TARGET_LINK_LIBRARIES( test_pic_accumulator
#     sll_pic_utilities
#     sll_pic_initializers
#     sll_pic_accumulators
#     sll_pic_particles
#     sll_meshes
#     sll_assert )
#   ADD_TEST( NAME pic_accumulator COMMAND test_pic_accumulator )
#   SET_TESTS_PROPERTIES( pic_accumulator PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
# ENDIF(BUILD_TESTING)


