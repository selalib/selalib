
SET(ARGS " ")
SET(PROCS 4 CACHE STRING "Number of processors used for parallel tests" )
ADD_MPI_TEST(collective test_collective ${PROCS} ${ARGS})
SET_TESTS_PROPERTIES(collective PROPERTIES FAIL_REGULAR_EXPRESSION "NOT PASS")

ADD_MPI_TEST(remap_2d test_remap_2d ${PROCS} ${ARGS})
SET_TESTS_PROPERTIES(remap_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

ADD_MPI_TEST(remap_3d test_remap_3d ${PROCS} ${ARGS})
SET_TESTS_PROPERTIES(remap_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

IF(NOT STDF95)

   SET(PROCS 4)
   ADD_MPI_TEST( point_to_point_comms_1d test_p2p_comms_1d ${PROCS} ${ARGS})
   SET_TESTS_PROPERTIES( point_to_point_comms_1d PROPERTIES 
     PASS_REGULAR_EXPRESSION  "PASSED")

   SET(PROCS 4)
   ADD_MPI_TEST( point_to_point_comms_2d test_p2p_comms_2d ${PROCS} ${ARGS})
   SET_TESTS_PROPERTIES( point_to_point_comms_2d PROPERTIES 
     PASS_REGULAR_EXPRESSION "PASSED")

ENDIF(NOT STDF95)

IF(PROCESSOR_COUNT GREATER 1)

   ADD_MPI_TEST(remap_4d test_remap_4d ${PROCS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ADD_MPI_TEST(remap_5d test_remap_5d ${PROCS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_5d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

   ADD_MPI_TEST(remap_6d test_remap_6d ${PROCS} ${ARGS})
   SET_TESTS_PROPERTIES(remap_6d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


   IF(NOT STDF95)

      ADD_MPI_TEST(parallel_array_initializers test_parallel_array_initializer ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(parallel_array_initializers
	PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")


   ENDIF()

ENDIF(PROCESSOR_COUNT GREATER 1)


IF(HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL)
  
  IF(NOT STDF95)
    
    ADD_MPI_TEST(io_parallel test_io_parallel ${PROCS} ${ARGS})
    SET_TESTS_PROPERTIES(io_parallel PROPERTIES PASS_REGULAR_EXPRESSION 
      "PASSED")
    ######
    IF(FFT_LIB MATCHES "SLLFFT")
      
      ADD_MPI_TEST(poisson_3d_periodic_par 
        test_poisson_3d_periodic_par ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(poisson_3d_periodic_par 
        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
      
      ADD_MPI_TEST(poisson_per_cart_par_2d 
        test_poisson_2d_per_cart_par ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(poisson_per_cart_par_2d
        PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
      
      # SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/vpsim4d_input.txt)
      # ADD_MPI_TEST(vp4d_sim test_4d ${PROCS} ${ARGS})
      # SET_TESTS_PROPERTIES(vp4d_sim PROPERTIES 
      # PASS_REGULAR_EXPRESSION "PASSED")
      
      SET(ARGS ${CMAKE_BINARY_DIR}/vpsim4d_general_input.txt)
      ADD_MPI_TEST(vp4d_sim_general test_4d_vp_general ${PROCS} ${ARGS})

      SET_TESTS_PROPERTIES(vp4d_sim_general PROPERTIES 
	PASS_REGULAR_EXPRESSION "PASSED")
      ADD_MPI_TEST(qns2d_parallel test_qn_solver_2d_parallel ${PROCS} 
	${ARGS})
      SET_TESTS_PROPERTIES(qns2d_parallel PROPERTIES 
	PASS_REGULAR_EXPRESSION "PASSED")
      
    ENDIF()
    
    ######
    
    SET(ARGS ${CMAKE_BINARY_DIR}/vpsim4d_cartesian_input)
    ADD_MPI_TEST(sim4d_vp_cart test_4d_vp_cartesian ${PROCS} ${ARGS})
    SET_TESTS_PROPERTIES(sim4d_vp_cart PROPERTIES 
      PASS_REGULAR_EXPRESSION "PASSED")
    
    SET(ARGS ${CMAKE_BINARY_DIR}/vpsim2d_cartesian_input)
    ADD_MPI_TEST(sim2d_vp_cart test_2d_vp_cartesian ${PROCS} ${ARGS})
    SET_TESTS_PROPERTIES(sim2d_vp_cart PROPERTIES 
      PASS_REGULAR_EXPRESSION "PASSED")
    
    SET(ARGS ${CMAKE_BINARY_DIR}/dksim4d_polar_input.nml)
    ADD_MPI_TEST(sim4d_DK_polar test_4d_dk_polar ${PROCS} ${ARGS})
    SET_TESTS_PROPERTIES(sim4d_DK_polar PROPERTIES 
      PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 100)

    SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/dksim4d_polar_one_mu.nml)
    ADD_MPI_TEST(sim4d_DK_polar_one_mu test_4d_dk_polar_one_mu ${PROCS} ${ARGS})
    SET_TESTS_PROPERTIES(sim4d_DK_polar_one_mu PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 100)

    
    SET(ARGS ${CMAKE_BINARY_DIR}/sim4d_qns_general_input.txt)
    ADD_MPI_TEST(vp4d_sim_qns_general test_4d_qns_general ${PROCS} ${ARGS})
    SET_TESTS_PROPERTIES(vp4d_sim_qns_general PROPERTIES 
      PASS_REGULAR_EXPRESSION "PASSED")
   
    # SET(PROCS 1)
    # SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/sim4d_qns_mixed_input.txt)
    # ADD_MPI_TEST(vp4d_sim_mixed_qns_cartesian test_4d_mixed_qns_cartesian 
    # ${PROCS} ${ARGS})
    # SET_TESTS_PROPERTIES(vp4d_sim_mixed_qns_cartesian PROPERTIES 
    # PASS_REGULAR_EXPRESSION "PASSED")
    
    SET(ARGS ${CMAKE_BINARY_DIR}/sim4d_DK_hybrid_input.txt)
    ADD_MPI_TEST(sim4d_DK_hybrid test_4d_DK_hybrid ${PROCS} ${ARGS})
    SET_TESTS_PROPERTIES(sim4d_DK_hybrid PROPERTIES 
      PASS_REGULAR_EXPRESSION "PASSED")
    
    #SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/dksim4d_general_input.txt)
    ADD_MPI_TEST(dk4d_sim_cartesian test_4d_dk_cartesian ${PROCS} ${ARGS})
    SET_TESTS_PROPERTIES(dk4d_sim_cartesian PROPERTIES 
      PASS_REGULAR_EXPRESSION "PASSED")
    
    IF (FFT_LIB MATCHES "SLLFFT")
      SET(PROCS 4)
      #SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/dksim4d_general_input.txt)
      ADD_MPI_TEST(vp_4d_eulerian test_4d_vp_eulerian_cartesian_finite_volume 
	${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(vp_4d_eulerian PROPERTIES TIMEOUT 30 
	PASS_REGULAR_EXPRESSION "PASSED")
    ENDIF()
   
    IF(PROCESSOR_COUNT GREATER 1)
      
      SET(ARGS ${CMAKE_BINARY_DIR}/vpsim6d_input.txt)
      ADD_MPI_TEST(vp6d_sim test_6d ${PROCS} ${ARGS})
      SET_TESTS_PROPERTIES(vp6d_sim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
      
    ENDIF(PROCESSOR_COUNT GREATER 1)
    
  ENDIF(NOT STDF95)
  
  SET(ARGS " ")
  IF(MUDPACK_ENABLED AND Fortran_COMPILER STREQUAL "GFORTRAN")
    #ADD_MPI_TEST(multigrid_2d test_mgd2 ${PROCS} ${ARGS})
    #SET(multigrid_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
    ADD_MPI_TEST(multigrid_3d test_mgd3 ${PROCS} ${ARGS})
    SET(multigrid_3d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  ENDIF()
  
  SET(PROCS 1)
  ADD_MPI_TEST( visu_pic test_visu_pic ${PROCS} ${ARGS} )
  SET( visu_pic PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  SET(PROCS 1)
  SET(ARGS " ")
  ADD_MPI_TEST( pic_simulation_4d test_4d_vp_pic_cartesian ${PROCS} ${ARGS} )
  SET( pic_simulation_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")



ENDIF() # HDF5_PARALLEL_ENABLED AND HDF5_IS_PARALLEL

IF(PASTIX_FOUND AND PTSCOTCH_FOUND AND MURGE_FOUND AND SCOTCH_FOUND)
  SET(ARGS "1000 3") 
  ADD_MPI_TEST(pastix test_pastix ${PROCS} ${ARGS})
  SET(pastix PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDIF()
