ADD_DEFINITIONS(-DMUDPACK)

SET(POISSON_SRCS sll_m_poisson_1d_base.F90
                 sll_m_poisson_1d_periodic.F90
                 sll_m_poisson_1d_periodic_wrapper.F90
                 sll_m_poisson_1d_polar.F90
                 sll_m_poisson_1d_fourier.F90
                 sll_m_poisson_1d_fem.F90
                 sll_m_poisson_1d_fd.F90
                 sll_m_poisson_2d_polar.F90
                 sll_m_poisson_2d_polar_wrapper.F90
                 sll_m_poisson_2d_periodic_fftpack.F90
                 sll_m_poisson_2d_mudpack.F90
                 sll_m_poisson_2d_mudpack_curvilinear_solver_old.F90
                 sll_m_poisson_2d_base.F90
                 sll_m_poisson_2d_fft.F90
                 sll_m_poisson_2d_periodic_solver.F90
                 sll_m_poisson_2d_elliptic_solver.F90
                 sll_m_poisson_2d_curvilinear.F90
                 sll_m_poisson_2d_fem.F90
                 sll_m_poisson_3d_base.F90
                 sll_m_poisson_3d_periodic.F90
                 sll_m_fishpack.F90
                 sll_m_mudpack.F90
                 sll_m_mudpack_curvilinear.F90
                 sll_m_hex_poisson.F90
                 sll_m_tri_poisson.F90 )

IF(FFTW_ENABLED AND FFTW_FOUND)
   SET(POISSON_SRCS ${POISSON_SRCS} sll_m_poisson_2d_periodic_fftw.F90)
ENDIF()

ADD_LIBRARY( sll_poisson_solvers STATIC ${POISSON_SRCS} )
TARGET_LINK_LIBRARIES( sll_poisson_solvers
  sll_general_coordinate_elliptic_solvers 
  sll_meshes
  sll_memory 
  sll_assert 
  mudpack
  fishpack
  dfftpack 
  ${LAPACK_LIBRARIES})

# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
