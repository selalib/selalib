
ADD_LIBRARY( sll_general_coordinate_elliptic_solvers STATIC
     connectivity_qn.f90
     sll_general_coordinate_elliptic_solver.F90
     sll_knots.F90)

TARGET_LINK_LIBRARIES( sll_general_coordinate_elliptic_solvers 
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
                       sll_assert
                       sll_memory
                       sll_working_precision
                       deboor_splines)
   
ADD_EXECUTABLE(test_general_coordinates_elliptic_solver 
               init_functions.F90 unit_test.F90)

TARGET_LINK_LIBRARIES(test_general_coordinates_elliptic_solver
                       sll_general_coordinate_elliptic_solvers)

ADD_EXECUTABLE(test_gces_identity init_functions.F90
               unit_test_identity.F90)

TARGET_LINK_LIBRARIES(test_gces_identity sll_general_coordinate_elliptic_solvers)

ADD_EXECUTABLE(test_gces_colella_periodic init_functions.F90
               unit_test_colella_periodic.F90)

TARGET_LINK_LIBRARIES(test_gces_colella_periodic 
                      sll_general_coordinate_elliptic_solvers)

