SET(MAXWELL_SOLVERS_SRCS sll_m_maxwell_solvers_base.F90
                         sll_m_maxwell_2d_fdtd.F90
                         sll_m_maxwell_2d_diga.F90
                         sll_m_maxwell_1d_base.F90
                         sll_m_maxwell_1d_fem.F90)

IF(PGI)
   LIST(APPEND MAXWELL_SOLVERS_SRCS sll_m_maxwell_2d_pstd_pgi.F90
                                    sll_m_maxwell_3d_pstd_pgi.F90)
ELSE(PGI)
   LIST(APPEND MAXWELL_SOLVERS_SRCS sll_m_maxwell_2d_pstd.F90
                                    sll_m_maxwell_3d_pstd.F90)
ENDIF(PGI)

ADD_LIBRARY(sll_maxwell_solvers STATIC ${MAXWELL_SOLVERS_SRCS})

TARGET_LINK_LIBRARIES(sll_maxwell_solvers 
                      sll_fft
                      sll_fields 
                      sll_utilities)

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
