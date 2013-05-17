ADD_EXECUTABLE(VP_DG main.F90
                     mesh.F90
                     lob.f
                     gauss.f
                     mod_sparse.F90
                     vlasov_poisson_DG.F90
                   #  mod_octave_io.f90
                   #  mod_octave_io_sparse.f90
                     mod_umfpack.F90
                     umf4_f77wrapper.c)

TARGET_LINK_LIBRARIES(VP_DG sll_integration
                            sll_working_precision 
                            sll_constants 
                            ${SUITESPARSE_LIBRARIES} 
                            ${LAPACK_LIBRARIES})


IF(0)

   ADD_EXECUTABLE(test_DG test.F90
                          mesh.F90
                          gauss-lobatto.F90
                          lob.f
                          gauss.f
                          mod_sparse.F90
                        #  mod_octave_io.f90
                        #  mod_octave_io_sparse.f90
                 )

  ADD_DEPENDENCIES(test_DG sll_file_io)
  TARGET_LINK_LIBRARIES(test_DG sll_working_precision
	                             sll_constants)

ENDIF(0)
