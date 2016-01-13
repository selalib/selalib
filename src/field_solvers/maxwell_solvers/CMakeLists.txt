IF(FFTW_FOUND AND FFTW_ENABLED)

   ADD_LIBRARY(sll_maxwell_solvers STATIC sll_m_maxwell_solvers_base.F90
                                          sll_m_maxwell_2d_fdtd.F90
                                          sll_m_maxwell_2d_pstd.F90
                                          sll_m_maxwell_2d_diga.F90
                                          sll_m_maxwell_3d_pstd.F90
                                          sll_m_maxwell_1d_base.F90
                                          sll_m_maxwell_1d_fem.F90)

   TARGET_LINK_LIBRARIES(sll_maxwell_solvers 
                         sll_fft
                         sll_fields 
                         sll_utilities)

ELSE()

   SET(MAXWELL_SOLVERS_SRCS sll_m_maxwell_solvers_base.F90
                            sll_m_maxwell_1d_base.F90
                            sll_m_maxwell_1d_fem.F90
                            sll_m_maxwell_2d_fdtd.F90
                            sll_m_maxwell_2d_diga.F90)

   ADD_LIBRARY(sll_maxwell_solvers STATIC ${MAXWELL_SOLVERS_SRCS})

   TARGET_LINK_LIBRARIES(sll_maxwell_solvers sll_fields 
                                          sll_utilities )
ENDIF(FFTW_FOUND AND FFTW_ENABLED)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
