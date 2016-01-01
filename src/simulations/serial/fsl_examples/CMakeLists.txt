ADD_EXECUTABLE(fsl_advection_2d advection_2d.F90)
TARGET_LINK_LIBRARIES(fsl_advection_2d sll_interpolators )

#ADD_EXECUTABLE(vlasov_poisson_2d vlasov_poisson_2d.F90)
#TARGET_LINK_LIBRARIES(vlasov_poisson_2d sll_interpolators nufft )

ADD_EXECUTABLE(fsl_landau_2d landau_2d.F90)
TARGET_LINK_LIBRARIES(fsl_landau_2d sll_interpolators 
                                    sll_poisson_solvers
                                    sll_maxwell_solvers )

IF(FFTW_FOUND AND FFTW_ENABLED)
  ADD_EXECUTABLE(test_fftw3 test_fftw3.F90)
  TARGET_LINK_LIBRARIES(test_fftw3 sll_fft ${FFTW_LIBRARIES} )
ENDIF(FFTW_FOUND AND FFTW_ENABLED)
