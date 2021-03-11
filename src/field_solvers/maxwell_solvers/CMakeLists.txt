set(MAXWELL_SOLVERS_SRCS
    sll_m_maxwell_solvers_base.F90 sll_m_maxwell_2d_fdtd.F90
    sll_m_maxwell_2d_diga.F90 sll_m_maxwell_1d_base.F90
    sll_m_maxwell_1d_fem.F90)

list(APPEND MAXWELL_SOLVERS_SRCS sll_m_maxwell_2d_pstd.F90
     sll_m_maxwell_3d_pstd.F90)

add_library(sll_maxwell_solvers STATIC ${MAXWELL_SOLVERS_SRCS})

target_link_libraries(sll_maxwell_solvers sll_fft sll_fields sll_utilities)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
