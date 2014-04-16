IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)

ADD_EXECUTABLE(test_4d_vp_eulerian_cartesian_finite_volume 
               simulation_4d_vp_eulerian_cartesian_finite_volume.F90
               unit_test_4d_vp_eulerian_cartesian_finite_volume.F90)

TARGET_LINK_LIBRARIES(test_4d_vp_eulerian_cartesian_finite_volume 
                      selalib
                      ${HDF5_LIBRARIES})

IF(LOBALAP)
ADD_LIBRARY(sll_lobatto sll_lobatto_poisson.F90 lobalap.F90)

ADD_EXECUTABLE(test_lobalap unit_test_lobalap.F90 )
TARGET_LINK_LIBRARIES(test_lobalap sll_lobatto selalib ${HDF5_LIBRARIES})
ENDIF(LOBALAP)

ENDIF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)
