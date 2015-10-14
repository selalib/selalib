ADD_LIBRARY(sll_particle_pushers STATIC sll_operator_splitting_pic_2d2v_vp.F90 
  sll_operator_splitting_pic_1d2v_vm.F90
  sll_operator_splitting_cef_pic_1d2v_vm.F90)

TARGET_LINK_LIBRARIES( sll_particle_pushers
  sll_memory
  sll_assert
  sll_operator_splitting
  sll_pic_interface
  sll_kernel_smoothers
  sll_poisson_solvers
  sll_maxwell_solvers
  )


## Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF()
