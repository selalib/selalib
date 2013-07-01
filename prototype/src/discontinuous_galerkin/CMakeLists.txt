add_executable(VP_DG
	main.F90
  mesh.F90
  ../integration/sll_gausslobatto.F90
  ../integration/lob.f
  ../integration/gauss.f
  mod_sparse.F90
  sll_poisson_1d_periodic_dg.F90
  sll_timestep4dg.F90
  mod_octave_io.f90
  mod_octave_io_sparse.f90
  mod_umfpack.F90
  umf4_f77wrapper.c
	)
#add_dependencies(VP_DG sll_file_io )
target_link_libraries(VP_DG
  sll_working_precision
	sll_constants
	sll_utilities
  umfpack
  amd
  suitesparseconfig
  cholmod
  colamd
  ccolamd
  rt
  blas
  lapack
	)
add_dependencies(VP_DG sll_file_io )

if(0)
#this is a serie of test I used while developping
#this should disappear or be transformed into several unit_test
add_executable(test_DG
	test.F90
  mesh.F90
  ../integration/gauss-lobatto.F90
  lob.f
  gauss.f
  mod_sparse.F90
  mod_octave_io.f90
  mod_octave_io_sparse.f90
	)
add_dependencies(test_DG sll_file_io)
target_link_libraries(test_DG
  sll_working_precision
	sll_constants
	)
endif(0)
