IF (HDF5_IS_PARALLEL)
  SET( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF (HDF5_IS_PARALLEL)

ADD_EXECUTABLE( test_general_coordinates_elliptic_solver 
                m_init_functions.F90 test_general_coordinates_elliptic_solver.F90 )

TARGET_INCLUDE_DIRECTORIES(test_general_coordinates_elliptic_solver BEFORE PRIVATE ${HDF5_INCLUDE_DIR_FORTRAN})
TARGET_LINK_LIBRARIES( test_general_coordinates_elliptic_solver
                       sll_general_coordinate_elliptic_solvers )

ADD_TEST( NAME general_coordinate_elliptic_solver 
          COMMAND test_general_coordinates_elliptic_solver )

SET_TESTS_PROPERTIES( general_coordinate_elliptic_solver 
                      PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 2000 )

ADD_EXECUTABLE( test_lobalap_analytic test_lobalap_analytic.F90 )
TARGET_INCLUDE_DIRECTORIES(test_lobalap_analytic BEFORE PRIVATE ${HDF5_INCLUDE_DIR_FORTRAN})
TARGET_LINK_LIBRARIES( test_lobalap_analytic sll_general_coordinate_elliptic_solvers )

IF(NOT PGI)

  ADD_EXECUTABLE( test_lobalap_discrete test_lobalap_discrete.F90 )
  TARGET_INCLUDE_DIRECTORIES(test_lobalap_discrete BEFORE PRIVATE ${HDF5_INCLUDE_DIR_FORTRAN})
  TARGET_LINK_LIBRARIES( test_lobalap_discrete sll_general_coordinate_elliptic_solvers )

  ADD_EXECUTABLE( test_gces_full_periodic_prototype test_full_periodic_prototype.F90 )
  TARGET_INCLUDE_DIRECTORIES(test_gces_full_periodic_prototype BEFORE PRIVATE ${HDF5_INCLUDE_DIR_FORTRAN})
  TARGET_LINK_LIBRARIES( test_gces_full_periodic_prototype sll_general_coordinate_elliptic_solvers )

ENDIF(NOT PGI)
