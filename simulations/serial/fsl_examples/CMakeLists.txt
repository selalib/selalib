ADD_EXECUTABLE(fsl_advection_2d advection_2d.F90)
TARGET_LINK_LIBRARIES(fsl_advection_2d sll_interpolators )

ADD_EXECUTABLE(fsl_vlasov_poisson_2d deposit_cubic_splines.F90
                                     vlasov_poisson_2d.F90)

TARGET_LINK_LIBRARIES(fsl_vlasov_poisson_2d sll_interpolators sll_fft )

ADD_EXECUTABLE(fsl_landau_2d landau_2d.F90)
TARGET_LINK_LIBRARIES(fsl_landau_2d sll_interpolators sll_poisson_solvers)
