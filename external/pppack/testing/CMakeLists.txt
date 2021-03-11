add_executable(test_pppack test_pppack.F90)
target_link_libraries(test_pppack pppack)
add_executable(test_woodbury_solver test_woodbury_solver.F90)
target_link_libraries(test_woodbury_solver pppack)
add_executable(test_schur_complement_solver test_schur_complement_solver.F90)
target_link_libraries(test_schur_complement_solver pppack)

# FIND_PACKAGE(GSL) IF(GSL_FOUND) ADD_EXECUTABLE( test_gsl_bsplines
# test_gsl_bsplines.c ) TARGET_LINK_LIBRARIES( test_gsl_bsplines
# ${GSL_LIBRARIES}) ENDIF(GSL_FOUND)
