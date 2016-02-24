ADD_LIBRARY(sll_pic_time_integration STATIC 
  sll_m_operator_splitting_pic_vp_2d2v.F90 
  sll_m_hamiltonian_splitting_pic_vm_1d2v.F90)

TARGET_LINK_LIBRARIES( sll_pic_time_integration
  sll_memory
  sll_assert
  sll_splitting_methods
  sll_pic_interface
  sll_kernel_smoothers
  sll_poisson_solvers
  sll_maxwell_solvers
  sll_collective
  )


## Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF()
