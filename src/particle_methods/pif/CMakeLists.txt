add_library(sll_pif STATIC sll_m_pif_fieldsolver.F90
                           sll_m_wedge_product_generaldim.F90)
target_link_libraries(sll_pif sll_memory sll_assert sll_collective
                      sll_utilities)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
