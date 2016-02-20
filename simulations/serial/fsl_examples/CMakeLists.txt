ADD_EXECUTABLE(fsl_advection_2d advection_2d.F90)
TARGET_LINK_LIBRARIES(fsl_advection_2d sll_interpolators )

IF(FFTW_FOUND)
  ADD_EXECUTABLE(fsl_vlasov_poisson_2d vlasov_poisson_2d.F90)
  TARGET_LINK_LIBRARIES(fsl_vlasov_poisson_2d sll_interpolators nufft ${FFTW_LIBRARIES} )
ENDIF(FFTW_FOUND)

ADD_EXECUTABLE(fsl_landau_2d landau_2d.F90)
TARGET_LINK_LIBRARIES(fsl_landau_2d sll_interpolators sll_poisson_solvers)
