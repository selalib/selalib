IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)

ADD_EXECUTABLE(test_4d_vp_eulerian_cartesian_finite_volume 
               simulation_4d_vp_eulerian_cartesian_finite_volume.F90
               unit_test_4d_vp_eulerian_cartesian_finite_volume.F90)

TARGET_LINK_LIBRARIES(test_4d_vp_eulerian_cartesian_finite_volume 
                      selalib
                      ${HDF5_LIBRARIES})

ADD_EXECUTABLE(lobalap lobalap.f90)
TARGET_LINK_LIBRARIES(lobalap sll_utilities) 

ENDIF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)
