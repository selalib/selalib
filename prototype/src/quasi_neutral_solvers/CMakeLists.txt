set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

#add_library(sll_quasi_neutral_solver STATIC sll_quasi_neutral_solver.F90 bsplvb.F90 bsplvd.F90)
#add_executable(test_quasi_neutral_solver unit_test.F90)
#target_link_libraries(test_quasi_neutral_solver sll_assertion sll_constants 
#sll_diagnostics sll_memory sll_mesh_types sll_utilities sll_remap sll_collective sll_diagnostics sll_timer sll_fft #sll_quasi_neutral_solver)


add_library(sll_qns2d_seq STATIC sll_qns2d_with_finite_diff_seq.F90
                                 sll_qns2d_angular_spectral_method_seq.F90)

add_dependencies(sll_qns2d_seq  
                 sll_fft 
                 sll_tridiagonal
                 sll_memory 
                 sll_constants 
                 sll_utilities 
                 sll_working_precision
                 sll_timer 
                 sll_assertion )

add_library(sll_qns2d_par STATIC sll_qns2d_with_finite_diff_par.F90
                                 sll_qns2d_angular_spectral_method_par.F90)


set(QNS_LIBS     sll_fft 
                 sll_numerical_utilities
                 sll_memory 
                 sll_constants 
                 sll_utilities 
                 sll_working_precision
                 sll_timer 
                 sll_assertion )

add_dependencies(sll_qns2d_par  
                 sll_remap   
                 sll_collective 
                 ${QNS_LIBS})

# 2D QNS with finite defferences method in sequential
add_executable( test_qns2d_with_finite_diff_seq
                test_sll_qns2d_with_finite_diff_seq.F90)

target_link_libraries( test_qns2d_with_finite_diff_seq sll_qns2d_seq
                                                       ${QNS_LIBS} )

# 2D with angular spectral method in sequential
add_executable( test_qns2d_angular_spectral_method_seq 
                test_sll_qns2d_angular_spectral_method_seq.F90)

target_link_libraries(test_qns2d_angular_spectral_method_seq sll_qns2d_seq 
                                                             ${QNS_LIBS} )


# 2D QNS with finite defferences method in parallel
add_executable( test_qns2d_with_finite_diff_par 
                test_sll_qns2d_with_finite_diff_par.F90)

target_link_libraries( test_qns2d_with_finite_diff_par sll_qns2d_par
                                                       sll_remap
                                                       sll_collective 
                                                       ${QNS_LIBS} )

# 2D QNS with angular spectral method in parallel
add_executable( test_qns2d_angular_spectral_method_par 
                test_sll_qns2d_angular_spectral_method_par.F90)

target_link_libraries( test_qns2d_angular_spectral_method_par sll_qns2d_par
                                                              sll_remap
                                                              sll_collective
                                                              ${QNS_LIBS} )

