set(POISSON_SRCS
    sll_m_poisson_1d_base.F90
    sll_m_poisson_1d_periodic.F90
    sll_m_poisson_1d_hmf.F90
    sll_m_poisson_1d_polar.F90
    sll_m_poisson_2d_polar.F90
    sll_m_poisson_2d_base.F90
    sll_m_poisson_2d_periodic.F90
    sll_m_poisson_2d_fem.F90
    sll_m_poisson_2d_fem_fft.F90
    sll_m_hex_poisson.F90
    sll_m_poisson_2d_tri.F90
    sll_m_poisson_2d_fem_sps_weak_form.F90
    sll_m_poisson_2d_fem_sps_dense_assembler.F90
    sll_m_poisson_2d_fem_sps_dense_projector.F90
    sll_m_poisson_2d_fem_sps_dense.F90
    # sll_m_poisson_2d_fem_sps_stencil_assembler.F90
    # sll_m_poisson_2d_fem_sps_stencil_projector.F90
    # sll_m_poisson_2d_fem_sps_stencil.F90
    sll_m_poisson_2d_fem_sps_stencil_new_assembler.F90
    sll_m_poisson_2d_fem_sps_stencil_new.F90
    sll_m_poisson_3d_base.F90
    sll_m_poisson_3d_fem_fft.F90    
    sll_m_preconditioner_poisson_fft.F90)

if(FFT_LIB MATCHES FFTW)
  list(APPEND POISSON_SRCS sll_m_poisson_3d_periodic.F90)
endif()

if(BUILD_MUDPACK)
  add_definitions(-DMUDPACK)
  set(POISSON_SRCS
      ${POISSON_SRCS} sll_m_poisson_2d_mudpack_curvilinear.F90
      sll_m_poisson_2d_mudpack_curvilinear_old.F90 sll_m_mudpack.F90
      sll_m_mudpack_curvilinear.F90)
endif(BUILD_MUDPACK)

add_library(sll_poisson_solvers STATIC ${POISSON_SRCS})

target_link_libraries(
  sll_poisson_solvers
  sll_assert
  sll_errors
  sll_memory
  sll_working_precision
  sll_boundary_condition_descriptors
  sll_fft
  sll_general_coordinate_elliptic_solvers
  sll_meshes
  sll_splines
  sll_fem_spline
  sll_splines_smooth
  ${LAPACK_LIBRARIES})

if(BUILD_MUDPACK)
  target_link_libraries(sll_poisson_solvers mudpack)
endif(BUILD_MUDPACK)

# Ctests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
