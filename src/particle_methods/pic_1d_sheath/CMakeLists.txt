ADD_LIBRARY( sll_pic_1d_sheath_field_solver STATIC pic_1d_field_solver.F90
  particle_description.F90 )

TARGET_LINK_LIBRARIES( sll_pic_1d_sheath_field_solver
  sll_constants
  sll_working_precision
  sll_memory
  sll_meshes
  sll_boundary_condition_descriptors
  sll_fft
  sll_splines
  sll_integration
  sll_utilities
  sll_collective
  sll_poisson_solvers )


ADD_LIBRARY( sll_pic_1d_sheath STATIC particle_loading.F90
  pic_1d.F90
  unit_test.F90
  pic_postprocessing.F90
  particle_description.F90 )

TARGET_LINK_LIBRARIES( sll_pic_1d_sheath
  sll_pic_1d_sheath_field_solver
  sll_constants
  sll_working_precision
  sll_memory
  sll_meshes
  sll_boundary_condition_descriptors
  sll_splines
  sll_integration
  sll_utilities
  sll_interpolators
  sll_collective
  sll_prob
  sll_sobol
  sll_visu_pic
  sll_timer )


ADD_EXECUTABLE( test_pic_1d_sheath unit_test.F90 )
TARGET_LINK_LIBRARIES( test_pic_1d_sheath sll_collective sll_pic_1d_sheath )

