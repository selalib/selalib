ADD_EXECUTABLE(bsl_1d_cubic_uniform_periodic test_1d_cubic_uniform_periodic.F90)
TARGET_LINK_LIBRARIES(bsl_1d_cubic_uniform_periodic sll_interpolators)

ADD_EXECUTABLE(bsl_1d_cubic_uniform_compact test_1d_cubic_uniform_compact.F90)
TARGET_LINK_LIBRARIES(bsl_1d_cubic_uniform_compact sll_interpolators)

ADD_EXECUTABLE(bsl_1d_cubic_nonuniform_periodic test_1d_cubic_nonuniform_periodic.F90)
TARGET_LINK_LIBRARIES(bsl_1d_cubic_nonuniform_periodic sll_interpolators)

ADD_EXECUTABLE(bsl_1d_cubic_nonuniform_compact test_1d_cubic_nonuniform_compact.F90)
TARGET_LINK_LIBRARIES(bsl_1d_cubic_nonuniform_compact sll_interpolators)

ADD_EXECUTABLE(bsl_1d_quintic_uniform_compact test_1d_quintic_uniform_compact.F90)
TARGET_LINK_LIBRARIES(bsl_1d_quintic_uniform_compact sll_interpolators)

ADD_EXECUTABLE(bsl_1d_quintic_nonuniform_compact test_1d_quintic_nonuniform_compact.F90)
TARGET_LINK_LIBRARIES(bsl_1d_quintic_nonuniform_compact sll_interpolators)

ADD_EXECUTABLE(bsl_2d_omp unit_test_2d_omp.F90)
ADD_DEPENDENCIES(bsl_2d_omp sll_interpolators)
TARGET_LINK_LIBRARIES(bsl_2d_omp sll_interpolators)

ADD_EXECUTABLE(bsl_landau_1d landau_1d.F90)
TARGET_LINK_LIBRARIES(bsl_landau_1d selalib dfftpack)

ADD_EXECUTABLE(bsl_landau_2d landau_2d.F90)
TARGET_LINK_LIBRARIES(bsl_landau_2d sll_interpolators 
                                    sll_poisson_solvers
                                    sll_maxwell_solvers )
#                                    mudpack  )
