
ADD_LIBRARY( sll_general_coordinate_elliptic_solvers STATIC
     connectivity_qn.f90
     sll_general_coordinate_elliptic_solver.F90
     sll_knots.F90)

ADD_LIBRARY( finite_elements_solver STATIC
     connectivity_qn.f90
     finite_elements_solver.F90
     sll_knots.F90)

TARGET_LINK_LIBRARIES( sll_general_coordinate_elliptic_solvers
		       finite_elements_solver 
		       sll_logical_meshes
                       sll_sparse_matrix 
                       sll_scalar_field_2d_alternative 
                       sll_interpolators
                       sll_coordinate_transformations
                       sll_logical_meshes
                       sll_boundary_condition_descriptors
                       sll_file_io
                       sll_integration
                       sll_constants
                       sll_utilities
                       sll_timer
                       sll_assertion
                       sll_memory
                       sll_working_precision
                       deboor_splines)
   
ADD_EXECUTABLE(test_general_coordinates_elliptic_solver unit_test.F90)

TARGET_LINK_LIBRARIES( test_general_coordinates_elliptic_solver
                       sll_general_coordinate_elliptic_solvers)

ADD_EXECUTABLE(test_finite_elements_solver unit_test_finite_elements.F90)

TARGET_LINK_LIBRARIES( test_finite_elements_solver
                       sll_general_coordinate_elliptic_solvers)
