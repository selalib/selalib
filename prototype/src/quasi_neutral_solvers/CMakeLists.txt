add_library(sll_quasi_neutral_solver STATIC sll_quasi_neutral_solver.F90 bsplvb.F90 bsplvd.F90)
add_executable(test_quasi_neutral_solver unit_test.F90)
target_link_libraries(test_quasi_neutral_solver sll_assertion sll_constants 
 sll_diagnostics sll_memory sll_mesh_types sll_utilities sll_remap sll_collective sll_diagnostics sll_timer fftpack_module sll_quasi_neutral_solver)


# Aliou Codes

add_library(sll_qns2d STATIC sll_qns2d_angular_spect_method_seq.F90 
                             sll_qns2d_angular_spect_method_par.F90
                                sll_qns_2d_with_finite_diff_seq.F90
                               sll_qns_2d_with_finite_diff.F90_par)


add_executable(test_qns2d test_qns2d.F90)

target_link_libraries(test_qns2d sll_assertion sll_constants 
                   sll_diagnostics sll_memory sll_mesh_types 
                      sll_utilities sll_remap sll_collective 
                                   sll_diagnostics sll_timer 
                    fftpack_module sll_qns2d)
