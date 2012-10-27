#add_library(sll_quasi_neutral_solver STATIC sll_quasi_neutral_solver.F90 bsplvb.F90 bsplvd.F90)
#add_executable(test_quasi_neutral_solver unit_test.F90)
#target_link_libraries(test_quasi_neutral_solver sll_assertion sll_constants 
#sll_diagnostics sll_memory sll_mesh_types sll_utilities sll_remap sll_collective sll_diagnostics sll_timer sll_fft #sll_quasi_neutral_solver)


add_library( sll_qns2d STATIC sll_qns2d_with_finite_diff_seq.F90
                              sll_qns2d_with_finite_diff_par.F90
                              sll_qns2d_angular_spectral_method_seq.F90
                              sll_qns2d_angular_spectral_method_par.F90
           )

set(QNS_LIBS sll_qns2d sll_assertion sll_fft 
             sll_memory sll_constants
             sll_utilities sll_remap   
             sll_collective sll_working_precision
             sll_timer sll_tridiagonal
    ${QNS_LIBS})

add_dependencies( ${QNS_LIBS} )

# 2D QNS with finite defferences method in sequential
add_executable( test_qns2d_with_finite_diff_seq
                test_sll_qns2d_with_finite_diff_seq.F90)
target_link_libraries( test_qns2d_with_finite_diff_seq ${QNS_LIBS} )

# 2D QNS with finite defferences method in parallel
add_executable( test_qns2d_with_finite_diff_par 
                test_sll_qns2d_with_finite_diff_par.F90)
target_link_libraries( test_qns2d_with_finite_diff_par ${QNS_LIBS} )

# 2D with angular spectral method in sequential
add_executable( test_qns2d_angular_spectral_method_seq 
                test_sll_qns2d_angular_spectral_method_seq.F90)
target_link_libraries(test_qns2d_angular_spectral_method_seq ${QNS_LIBS} )

# 2D QNS with angular spectral method in parallel
add_executable( test_qns2d_angular_spectral_method_par 
               test_sll_qns2d_angular_spectral_method_par.F90)
target_link_libraries( test_qns2d_angular_spectral_method_par ${QNS_LIBS} )

