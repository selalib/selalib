IF(FFTW_LIBRARY AND HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)

   ADD_EXECUTABLE(test_4d_vp_eulerian_cartesian_finite_volume 
                  simulation_4d_vp_eulerian_cartesian_finite_volume.F90
                  unit_test_4d_vp_eulerian_cartesian_finite_volume.F90)

   TARGET_LINK_LIBRARIES(test_4d_vp_eulerian_cartesian_finite_volume 
                         sll_simulations
                         sll_parallel_array_initializers
                         sll_coordinate_transformations
                         sll_poisson_parallel
                         sll_point_to_point_communications 
                         sll_logical_meshes
                         sll_dg_fields
                         sll_mesh_calculus
                         sll_integration
                         sll_file_io_parallel
                         sll_timer
                         sll_collective
                         ${FFTW_LIBRARY}
                         ${HDF5_LIBRARIES})

   ADD_LIBRARY(sll_lobatto sll_lobatto_poisson.F90 
                           map_function_module.F90
                           lobalap.F90)

   ADD_EXECUTABLE(test_lobalap unit_test_lobalap.F90 )

   TARGET_LINK_LIBRARIES(test_lobalap sll_lobatto 
                                      sll_dg_fields
                                      sll_coordinate_transformations
                                      sll_file_io)

   ADD_EXECUTABLE(test_lobalap_discrete unit_test_lobalap_discrete.F90 )

   TARGET_LINK_LIBRARIES(test_lobalap_discrete sll_lobatto 
                                               sll_dg_fields
                                               sll_file_io)

#    IF (FFT_LIB MATCHES "SLLFFT")
#      SET(PROCS 4)
#      #SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/dksim4d_general_input.txt)
#      ADD_MPI_TEST(vp_4d_eulerian test_4d_vp_eulerian_cartesian_finite_volume 
#	${PROCS} ${ARGS})
#      SET_TESTS_PROPERTIES(vp_4d_eulerian PROPERTIES TIMEOUT 30 
#	PASS_REGULAR_EXPRESSION "PASSED")
#    ENDIF()

ENDIF()
