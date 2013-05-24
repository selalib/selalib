ADD_LIBRARY( sll_general_coordinate_quasi_neutral_solvers STATIC
  sll_general_coordinate_qn_solver.F90
  coocsr.f90
  sll_knots.F90

  sparse_matrix_def.F90
  quicksort.f90
  sparse_matrix_LM.F90 
  sparse_matrix_basic.F90
  connectivities_def.F90
 
  sparse_matrix_composed.F90
  connectivity.F90

  sparse_matrix_operations.F90 
  sparse_matrix_con.F90
  sparse_matrix.F90
  GC.f90
 tracelog.F90
  precision.f90
  connectivity_qn.f90
  profilmat.f90 )

ADD_DEPENDENCIES( sll_general_coordinate_quasi_neutral_solvers
  sll_scalar_field_2d
  sll_integration
  sll_memory
  sll_boundary_condition_descriptors
  sll_working_precision
  deboor_splines
  )
