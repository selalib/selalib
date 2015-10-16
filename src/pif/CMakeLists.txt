ADD_LIBRARY( sll_pif_fieldsolver STATIC sll_m_pif_fieldsolver.F90 )
TARGET_LINK_LIBRARIES( sll_pif_fieldsolver
  sll_memory
  sll_assert
  sll_working_precision
  sll_collective
  sll_constants
  sll_utilities)

ADD_LIBRARY( sll_wedge_product_generaldim STATIC sll_m_wedge_product_generaldim.F90 )
TARGET_LINK_LIBRARIES( sll_wedge_product_generaldim
  sll_memory
  sll_assert
  sll_working_precision)
  
 
					      

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
