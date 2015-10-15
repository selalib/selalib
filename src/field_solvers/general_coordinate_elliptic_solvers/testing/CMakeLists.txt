ADD_EXECUTABLE( test_general_coordinates_elliptic_solver 
                init_functions.F90 test_general_coordinates_elliptic_solver.F90 )
TARGET_LINK_LIBRARIES( test_general_coordinates_elliptic_solver
                       sll_general_coordinate_elliptic_solvers )
ADD_TEST( NAME general_coordinate_elliptic_solver 
          COMMAND test_general_coordinates_elliptic_solver )
SET_TESTS_PROPERTIES( general_coordinate_elliptic_solver 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 2000 )

ADD_EXECUTABLE( test_lobalap_analytic test_lobalap_analytic.F90 )
TARGET_LINK_LIBRARIES( test_lobalap_analytic sll_general_coordinate_elliptic_solvers )

ADD_EXECUTABLE( test_lobalap_discrete test_lobalap_discrete.F90 )
TARGET_LINK_LIBRARIES( test_lobalap_discrete sll_general_coordinate_elliptic_solvers )

ADD_EXECUTABLE( test_gces_dirichlet_homogeneous test_dirichlet_homogeneous.F90 )
TARGET_LINK_LIBRARIES( test_gces_dirichlet_homogeneous sll_general_coordinate_elliptic_solvers )

ADD_TEST( NAME gces_dirichlet_homogeneous COMMAND test_gces_dirichlet_homogeneous )
SET_TESTS_PROPERTIES( gces_dirichlet_homogeneous PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 100 )

ADD_EXECUTABLE( test_gces_full_periodic test_full_periodic.F90 )
TARGET_LINK_LIBRARIES( test_gces_full_periodic sll_general_coordinate_elliptic_solvers )
#ADD_TEST( NAME gces_full_periodic COMMAND test_gces_full_periodic )
#SET_TESTS_PROPERTIES( gces_full_periodic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 100 )

SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ADD_EXECUTABLE( test_gces_dirichlet test_dirichlet.F90 )
TARGET_LINK_LIBRARIES( test_gces_dirichlet sll_general_coordinate_elliptic_solvers )

ADD_EXECUTABLE( test_gces_full_periodic_prototype test_full_periodic_prototype.F90 )
TARGET_LINK_LIBRARIES( test_gces_full_periodic_prototype sll_general_coordinate_elliptic_solvers )
