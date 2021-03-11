add_executable(pif_fieldsolver test_pif_fieldsolver.F90)
target_link_libraries(pif_fieldsolver sll_pif sll_timer)

add_executable(test_wedge_product_generaldim test_wedge_product_generaldim.F90)
target_link_libraries(test_wedge_product_generaldim sll_pif)

add_test(NAME pif_wedge_product_generaldim
         COMMAND test_wedge_product_generaldim)
set_tests_properties(pif_wedge_product_generaldim
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
