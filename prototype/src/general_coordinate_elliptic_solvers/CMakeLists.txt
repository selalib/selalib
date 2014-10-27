ADD_LIBRARY( sll_general_coordinate_elliptic_solvers STATIC
             connectivity_qn.f90
             sll_general_coordinate_elliptic_solver.F90
             sll_knots.F90
             sll_lobatto_poisson.F90 
             map_function_module.F90
             lobalap.F90)

TARGET_LINK_LIBRARIES( sll_general_coordinate_elliptic_solvers 
                       sll_sparse_matrix 
                       sll_scalar_field_2d_alternative 
                       sll_interpolators
                       sll_coordinate_transformations
                       sll_logical_meshes
                       sll_boundary_condition_descriptors
                       sll_dg_fields
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
               example_identity.F90)

TARGET_LINK_LIBRARIES(test_gces_identity sll_general_coordinate_elliptic_solvers)

ADD_EXECUTABLE(test_gces_colella_periodic init_functions.F90
               example_colella_periodic.F90)

TARGET_LINK_LIBRARIES(test_gces_colella_periodic 
                      sll_general_coordinate_elliptic_solvers)

ADD_TEST(NAME general_coordinate_elliptic_solver 
         COMMAND test_general_coordinates_elliptic_solver)

SET_TESTS_PROPERTIES(general_coordinate_elliptic_solver 
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_EXECUTABLE(test_lobalap_analytic unit_test_lobalap_analytic.F90 )

TARGET_LINK_LIBRARIES(test_lobalap_analytic sll_general_coordinate_elliptic_solvers )

ADD_EXECUTABLE(test_lobalap_discrete unit_test_lobalap_discrete.F90 )

TARGET_LINK_LIBRARIES(test_lobalap_discrete sll_general_coordinate_elliptic_solvers )
