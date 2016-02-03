ADD_DEFINITIONS(-DMUDPACK)

SET(POISSON_SRCS sll_m_poisson_1d_base.F90
                 sll_m_poisson_1d_periodic.F90
                 sll_m_poisson_1d_polar.F90
                 sll_m_poisson_1d_fourier.F90
                 sll_m_poisson_1d_fem.F90
                 sll_m_poisson_1d_fd.F90
                 sll_m_poisson_2d_polar.F90
                 sll_m_poisson_2d_periodic_fftpack.F90
                 sll_m_poisson_2d_base.F90
                 sll_m_poisson_2d_fft.F90
                 sll_m_poisson_2d_periodic.F90
                 sll_m_poisson_2d_elliptic_solver.F90
                 sll_m_poisson_2d_fem.F90
                 sll_m_poisson_3d_base.F90
                 sll_m_poisson_3d_periodic.F90
                 sll_m_hex_poisson.F90
                 sll_m_poisson_2d_tri.F90 )

IF(BUILD_MUDPACK)
  SET(POISSON_SRCS ${POISSON_SRCS}
                 sll_m_poisson_2d_mudpack_curvilinear.F90
                 sll_m_poisson_2d_mudpack_curvilinear_old.F90
                 sll_m_mudpack.F90
                 sll_m_mudpack_curvilinear.F90)
ENDIF(BUILD_MUDPACK)

IF(BUILD_FISHPACK)
  SET(POISSON_SRCS ${POISSON_SRCS} sll_m_fishpack.F90)
ENDIF(BUILD_FISHPACK)

IF(FFTW_ENABLED AND FFTW_FOUND)
   SET(POISSON_SRCS ${POISSON_SRCS} sll_m_poisson_2d_periodic_fftw.F90)
ENDIF()

ADD_LIBRARY( sll_poisson_solvers STATIC ${POISSON_SRCS} )
TARGET_LINK_LIBRARIES( sll_poisson_solvers
  sll_general_coordinate_elliptic_solvers 
  sll_meshes
  sll_memory 
  sll_assert 
  dfftpack 
  ${LAPACK_LIBRARIES})

IF(BUILD_FISHPACK)
  TARGET_LINK_LIBRARIES( sll_poisson_solvers fishpack)
ENDIF(BUILD_FISHPACK)

IF(BUILD_MUDPACK)
  TARGET_LINK_LIBRARIES( sll_poisson_solvers mudpack)
ENDIF(BUILD_MUDPACK)

IF(FFTW_ENABLED AND FFTW_FOUND)
   SET(POISSON_SRCS ${POISSON_SRCS} sll_m_poisson_2d_periodic_fftw.F90)
ENDIF()

# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
