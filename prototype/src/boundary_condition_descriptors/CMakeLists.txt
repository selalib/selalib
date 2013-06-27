add_library(sll_boundary_condition_descriptors STATIC 
  sll_boundary_condition_descriptors.F90)

add_dependencies(sll_boundary_condition_descriptors
  sll_working_precision
  )
