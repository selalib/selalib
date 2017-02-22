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

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
