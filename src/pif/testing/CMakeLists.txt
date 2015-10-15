ADD_EXECUTABLE(pif_fieldsolver_test  pif_fieldsolver_test.F90)
TARGET_LINK_LIBRARIES(pif_fieldsolver_test sll_pif_fieldsolver
					      sll_timer)
