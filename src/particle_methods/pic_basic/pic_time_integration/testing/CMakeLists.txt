add_executable(test_operator_splitting_pic_vp_2d2v
               test_operator_splitting_pic_vp_2d2v.F90)
target_link_libraries(test_operator_splitting_pic_vp_2d2v
                      sll_pic_time_integration)

add_executable(test_hamiltonian_splitting_pic_vm_1d2v
               test_hamiltonian_splitting_pic_vm_1d2v.F90)
target_link_libraries(test_hamiltonian_splitting_pic_vm_1d2v
                      sll_pic_time_integration)

add_executable(test_hamiltonian_splitting_pic_vm_1d2v_boris
               test_hamiltonian_splitting_pic_vm_1d2v_boris.F90)
target_link_libraries(test_hamiltonian_splitting_pic_vm_1d2v_boris
                      sll_pic_time_integration)

if(FFT_LIB STREQUAL "FFTW")
  add_mpi_test(hamiltonian_splitting_pic_vm_1d2v
               test_hamiltonian_splitting_pic_vm_1d2v "1" "")
  set_tests_properties(hamiltonian_splitting_pic_vm_1d2v
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  add_mpi_test(hamiltonian_splitting_pic_vm_1d2v_boris
               test_hamiltonian_splitting_pic_vm_1d2v_boris "1" "")
  set_tests_properties(hamiltonian_splitting_pic_vm_1d2v_boris
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  add_mpi_test(operator_splitting_pic_vp_2d2v
               test_operator_splitting_pic_vp_2d2v "1" "")
  set_tests_properties(operator_splitting_pic_vp_2d2v
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endif()
