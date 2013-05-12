add_executable(VP_DG
	main.F90
  mesh.F90
  gauss-lobatto.F90
  lob.f
  gauss.f
  mod_sparse.F90
  vlasov_poisson_DG.F90
#  mod_octave_io.f90
#  mod_octave_io_sparse.f90
  mod_umfpack.F90
  umf4_f77wrapper.c
	)
add_dependencies(VP_DG sll_file_io )
target_link_libraries(VP_DG ${UMFPACK_LIBRARIES} 
                            sll_working_precision sll_constants sll_utilities)

if(0)
add_executable(test_DG
	test.F90
  mesh.F90
  gauss-lobatto.F90
  lob.f
  gauss.f
  mod_sparse.F90
#  mod_octave_io.f90
#  mod_octave_io_sparse.f90
	)
add_dependencies(test_DG sll_file_io)
target_link_libraries(test_DG
  sll_working_precision
	sll_constants
	)
endif(0)
