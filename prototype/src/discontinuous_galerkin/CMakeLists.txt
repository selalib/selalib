add_executable(VP_DG
	main.F90
  mesh.F90
  gauss-lobatto.F90
  lob.f
  gauss.f
	)
add_dependencies(VP_DG sll_file_io)
target_link_libraries(VP_DG
  sll_working_precision
	sll_constants
	)


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
