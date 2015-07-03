ADD_LIBRARY( sll_pic1d STATIC
  particle_description.F90
  particle_loading.F90
  pic_postprocessing.F90
  pic_1d_field_solver.F90 )

TARGET_LINK_LIBRARIES( sll_pic1d
  sll_constants
  sll_working_precision
  sll_memory
  sll_meshes
  sll_boundary_condition_descriptors
  sll_poisson_solvers 
  sll_fft
  sll_splines
  sll_integration
  sll_utilities
  sll_interpolators
  sll_collective
  sobol_lgpl
  probability_lgpl
  sll_visu_pic
  sll_timer )

