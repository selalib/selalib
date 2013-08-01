
SET(ARGS " ")
SET(PROCS 2)
ADD_MPI_TEST(collective test_collective ${PROCS} ${ARGS})
SET_TESTS_PROPERTIES(collective PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS")

SET(PROCS 8)
ADD_MPI_TEST(remap_2d test_remap_2d ${PROCS} ${ARGS})
SET_TESTS_PROPERTIES(remap_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

SET(PROCS 8)
ADD_MPI_TEST(remap_3d test_remap_3d ${PROCS} ${ARGS})
SET_TESTS_PROPERTIES(remap_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

IF(PROCESSOR_COUNT GREATER 1)

   SET(PROCS 16)
   ADD_MPI_TEST(remap_4d test_remap_4d ${PROCS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   SET(PROCS 8)
   ADD_MPI_TEST(remap_6d test_remap_6d ${PROCS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_6d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   IF(NOT STDF95)


      SET(PROCS 2)
      ADD_MPI_TEST(parallel_array_initializers test_parallel_array_initializer ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(parallel_array_initializers
	PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

      SET(PROCS 16)
      ADD_MPI_TEST(qns2d_parallel test_qn_solver_2d_parallel ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(qns2d_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ENDIF()

ENDIF(PROCESSOR_COUNT GREATER 1)


IF(HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL)
   
   IF(NOT STDF95)

      SET(PROCS 16)
      ADD_MPI_TEST(poisson_3d_periodic_par 
                   test_poisson_3d_periodic_par ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(poisson_3d_periodic_par 
                           PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

      SET(PROCS 4)
      ADD_MPI_TEST(io_parallel test_io_parallel ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(io_parallel PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

      SET(PROCS 1)
      ADD_MPI_TEST(poisson_per_cart_par_2d 
                   test_poisson_2d_per_cart_par ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(poisson_per_cart_par_2d
                           PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   
      SET(PROCS 8)
      SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/vpsim4d_input.txt)
      ADD_MPI_TEST(vp4d_sim test_4d ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(vp4d_sim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   
      SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/vpsim4d_general_input.txt)
      ADD_MPI_TEST(vp4d_sim_general test_4d_vp_general ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(vp4d_sim_general PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

      SET(PROCS 16)
      SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/sim4d_qns_general_input.txt)
      ADD_MPI_TEST(vp4d_sim_qns_general test_4d_qns_general ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(vp4d_sim_qns_general PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   
      #SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/dksim4d_general_input.txt)
      ADD_MPI_TEST(dk4d_sim_cartesian test_4d_dk_cartesian ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(dk4d_sim_cartesian PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   
   
      IF(PROCESSOR_COUNT GREATER 1)
   
      SET(PROCS 8)
      SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/vpsim6d_input.txt)
      ADD_MPI_TEST(vp6d_sim test_6d ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(vp6d_sim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
    
      ENDIF(PROCESSOR_COUNT GREATER 1)
   ENDIF(NOT STDF95)


   SET(ARGS " ")
   IF(MUDPACK_ENABLED AND Fortran_COMPILER STREQUAL "GFORTRAN")
      SET(PROCS 4)
      #ADD_MPI_TEST(multigrid_2d test_mgd2 ${PROCS} ${ARGS})
      #SET(multigrid_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
      ADD_MPI_TEST(multigrid_3d test_mgd3 ${PROCS} ${ARGS})
      SET(multigrid_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ENDIF()

ENDIF()

IF(PASTIX_FOUND AND PTSCOTCH_FOUND AND MURGE_FOUND AND SCOTCH_FOUND)
  SET(PROCS 4)
  SET(ARGS "1000 3") 
  ADD_MPI_TEST(pastix test_pastix ${PROCS} ${ARGS})
  SET(pastix PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()
