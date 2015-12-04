ADD_LIBRARY( sll_pif STATIC sll_m_pif_fieldsolver.F90
			     sll_m_wedge_product_generaldim.F90)
TARGET_LINK_LIBRARIES( sll_pif
  sll_memory
  sll_assert
  sll_collective
  sll_utilities)
					      

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
