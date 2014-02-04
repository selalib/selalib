IF(SUITESPARSE_FOUND)

   INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIRS})

   ADD_LIBRARY(sll_vp_dg mesh.F90
                      ../integration/sll_gausslobatto.F90
                      ../integration/lob.f
                      ../integration/gauss.f
                      mod_sparse.F90
                      sll_poisson_1d_periodic_dg.F90
                      sll_timestep4dg.F90
                      mod_umfpack.F90
                      umf4_f77wrapper.c)

   ADD_EXECUTABLE(test_vp_dg main.F90)

   TARGET_LINK_LIBRARIES(test_vp_dg sll_working_precision
                                    sll_vp_dg
                                    sll_constants
                                    sll_utilities
                                    sll_boundary_condition_descriptors
                                    ${SUITESPARSE_LIBRARIES}
                                    ${LAPACK_LIBRARIES}
                                    ${BLAS_LIBRARIES})

   IF(UNIX AND NOT APPLE)
      TARGET_LINK_LIBRARIES(test_vp_dg rt)
   ENDIF()

ENDIF()
