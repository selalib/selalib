ADD_LIBRARY( sll_pic_1d STATIC
  sll_m_particle_1d_description.F90
  sll_m_pic_1d_particle_loading.F90
  sll_m_pic_postprocessing.F90
  sll_m_pic_1d_field_solver.F90
  sll_m_pic_1d_distribution.F90 )

TARGET_LINK_LIBRARIES( sll_pic_1d
  sll_memory
  sll_meshes
  sll_poisson_solvers 
  sll_fft
  sll_splines
  sll_quadrature
  sll_utilities
  sll_interpolators
  sll_collective
  sll_prob
  sll_sobol
  sll_pic_visu
  sll_timer )
