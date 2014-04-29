IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)

   ADD_EXECUTABLE(test_4d_vp_eulerian_cartesian_finite_volume 
                  simulation_4d_vp_eulerian_cartesian_finite_volume.F90
                  unit_test_4d_vp_eulerian_cartesian_finite_volume.F90)

   TARGET_LINK_LIBRARIES(test_4d_vp_eulerian_cartesian_finite_volume 
                         sll_simulations
                         sll_parallel_array_initializers
                         sll_coordinate_transformations
                         sll_poisson_parallel
                         sll_point_to_point_communications 
                         sll_dg_fields
                         sll_mesh_calculus
                         sll_integration
                         sll_file_io_parallel
                         sll_timer
                         sll_collective
                         ${HDF5_LIBRARIES}
                         ${FFTW_LIBRARY})

   ADD_LIBRARY(sll_lobatto sll_lobatto_poisson.F90 
                           map_function_module.F90
                           lobalap.F90)

   ADD_EXECUTABLE(test_lobalap unit_test_lobalap.F90 )
   TARGET_LINK_LIBRARIES(test_lobalap sll_lobatto 
                                      sll_dg_fields
                                      sll_coordinate_transformations
                                      sll_file_io
                                      ${HDF5_LIBRARIES})

ENDIF()
