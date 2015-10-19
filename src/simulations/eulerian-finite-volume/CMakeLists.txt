IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_ENABLED AND BUILD_TESTING)

   ADD_EXECUTABLE(test_4d_vp_eulerian_cartesian_finite_volume
                  simulation_4d_vp_eulerian_cartesian_finite_volume.F90
                  unit_test_4d_vp_eulerian_cartesian_finite_volume.F90)

   TARGET_LINK_LIBRARIES(test_4d_vp_eulerian_cartesian_finite_volume
                         sll_sim_parallel
                         sll_parallel_array_utilities
                         sll_coordinate_transformations
                         sll_poisson_parallel
                         sll_point_to_point_communications
                         sll_meshes
                         sll_fields
                         sll_mesh_calculus
                         sll_integration
                         sll_file_io_parallel
                         sll_timer
                         sll_collective)

#    IF (FFT_LIB MATCHES "SLLFFT")
#      SET(PROCS 4)
#      #SET(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/simulation/dksim4d_general_input.txt)
#      ADD_MPI_TEST(vp_4d_eulerian test_4d_vp_eulerian_cartesian_finite_volume 
#	${PROCS} ${ARGS})
#      SET_TESTS_PROPERTIES(vp_4d_eulerian PROPERTIES TIMEOUT 30 
#	PASS_REGULAR_EXPRESSION "PASSED")
#    ENDIF()

ENDIF()
