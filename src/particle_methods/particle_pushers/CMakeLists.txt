ADD_LIBRARY(sll_particle_pushers STATIC sll_operator_splitting_pic_2x2v_vp.F90 
  sll_operator_splitting_pic_1d2v_vm.F90)

TARGET_LINK_LIBRARIES( sll_particle_pushers
  sll_memory
  sll_assert
  sll_operator_splitting
  sll_pic_interface
  sll_kernel_smoothers
  sll_poisson_solvers
  sll_maxwell_solvers
  )


# Tests
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_operator_splitting_pic_1d2v_vm unit_test_operator_splitting_pic_1d2v_vm.F90)
  TARGET_LINK_LIBRARIES(test_operator_splitting_pic_1d2v_vm sll_particle_pushers)
  ADD_TEST( NAME operator_splitting_pic_1d2v_vm COMMAND kernel_smoother_spline_1d)

ENDIF(BUILD_TESTING)