
SET(POISSON_SRCS sll_m_poisson_1d_base.F90
                 sll_m_poisson_1d_periodic.F90
                 sll_m_poisson_1d_polar.F90
                 sll_m_poisson_2d_polar.F90
                 sll_m_poisson_2d_base.F90
                 sll_m_poisson_2d_periodic.F90
                 sll_m_poisson_2d_elliptic_solver.F90
                 sll_m_poisson_2d_fem.F90
                 sll_m_poisson_3d_base.F90
                 sll_m_poisson_3d_periodic.F90
                 sll_m_hex_poisson.F90
                 sll_m_poisson_2d_tri.F90 )

IF(BUILD_MUDPACK)
  ADD_DEFINITIONS(-DMUDPACK)
  SET(POISSON_SRCS ${POISSON_SRCS}
                 sll_m_poisson_2d_mudpack_curvilinear.F90
                 sll_m_poisson_2d_mudpack_curvilinear_old.F90
                 sll_m_mudpack.F90
                 sll_m_mudpack_curvilinear.F90)
ENDIF(BUILD_MUDPACK)


ADD_LIBRARY( sll_poisson_solvers STATIC ${POISSON_SRCS} )
TARGET_LINK_LIBRARIES( sll_poisson_solvers
  sll_general_coordinate_elliptic_solvers 
  sll_meshes
  sll_fft 
  sll_memory 
  sll_assert 
  ${LAPACK_LIBRARIES})

IF(BUILD_MUDPACK)
  TARGET_LINK_LIBRARIES( sll_poisson_solvers mudpack)
ENDIF(BUILD_MUDPACK)

# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
