SET(BSL_LIB_DEPS sll_interpolators
                 sll_poisson_1d_periodic
                 sll_splines 
                 -ldfftpack
                 sll_memory 
	         sll_working_precision
                 sll_utilities 
                 sll_assertion)

ADD_EXECUTABLE(test_bsl_1d unit_test_1d.F90)
ADD_DEPENDENCIES(test_bsl_1d ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(test_bsl_1d ${BSL_LIB_DEPS})

ADD_EXECUTABLE(test_bsl_2d unit_test_2d.F90)
ADD_DEPENDENCIES(test_bsl_2d ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(test_bsl_2d ${BSL_LIB_DEPS})

ADD_EXECUTABLE(test_bsl_landau_2d landau_2d.F90)
ADD_DEPENDENCIES(test_bsl_landau_2d ${BSL_LIB_DEPS})
TARGET_LINK_LIBRARIES(test_bsl_landau_2d ${BSL_LIB_DEPS})

#ADD_LIBRARY(sll_bsl STATIC sll_bsl.F90)

#add_executable(test_landau_4d landau_4d.F90)
#target_link_libraries(test_landau_4d sll_bsl sll_poisson_2d_periodic sll_constants sll_splines sll_distribution_function sll_memory sll_mesh_types sll_diagnostics sll_utilities fftpack_module sll_assertion)
