add_executable(VP_DG
	main.F90
  mesh.F90
  gauss-lobatto.F90
	)
add_dependencies(VP_DG sll_file_io)
target_link_libraries(VP_DG
  sll_working_precision
	sll_memory
	sll_constants
	)
