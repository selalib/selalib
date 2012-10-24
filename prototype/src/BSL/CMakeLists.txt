ADD_LIBRARY(sll_bsl STATIC sll_bsl.F90)
ADD_EXECUTABLE(test_bsl unit_test.F90)
TARGET_LINK_LIBRARIES(test_bsl sll_bsl 
                               sll_fft 
                               sll_constants 
                               sll_splines 
                               sll_distribution_function 
                               sll_memory 
                               sll_utilities 
                               sll_assertion)

#add_executable(test_landau_4d landau_4d.F90)
#target_link_libraries(test_landau_4d sll_bsl sll_poisson_2d_periodic sll_constants sll_splines sll_distribution_function sll_memory sll_mesh_types sll_diagnostics sll_utilities fftpack_module sll_assertion)

#add_executable(test_landau_2d landau_2d.F90)
#target_link_libraries(test_landau_2d sll_bsl sll_poisson_1d_periodic sll_constants sll_splines sll_distribution_function sll_memory sll_mesh_types sll_diagnostics sll_utilities fftpack_module sll_assertion)
