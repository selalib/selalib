ADD_LIBRARY( deboor_splines STATIC
             sll_module_deboor_splines_2d.F90
             sll_module_deboor_splines_1d.F90 )

TARGET_LINK_LIBRARIES( deboor_splines 
                       sll_utilities
                       sll_memory
                       sll_working_precision
                       sll_assert
                       sll_boundary_condition_descriptors )
  
