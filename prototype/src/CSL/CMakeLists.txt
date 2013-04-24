IF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ENDIF()

add_library(sll_csl STATIC sll_csl.F90)

IF (NOT STDF95)
add_executable(test_csl unit_test.F90)
target_link_libraries(test_csl
  sll_csl 
  sll_distribution_function 
  sll_interpolators 
  sll_splines 
  sll_mapped_meshes 
  ode_solvers
  sll_scalar_field_2d 
  geometry_functions 
  sll_file_io 
  sll_memory
  sll_assertion 
  sll_working_precision 
  sll_utilities 
  sll_constants )
ENDIF(NOT STDF95)

add_library(sll_bgk_mesh STATIC sll_bgk_mesh.F90)

add_executable(vp_csl_bloc_unif vp_csl_bloc_unif.F90)
target_link_libraries(vp_csl_bloc_unif
	sll_bgk_mesh
	sll_csl
	sll_distribution_function
	sll_interpolators
	sll_constants
	sll_splines
	sll_file_io 
	sll_mapped_meshes
	sll_utilities
	ode_solvers
	sll_scalar_field_2d
	geometry_functions
	sll_assertion 
        sll_memory 
        sll_working_precision)

#add_library(sll_classic_csl STATIC sll_classic_csl.F90)

# Specific files which are not part of a library module?? What is going on here?
add_executable(vp_classic_csl bgk_classic_csl.F90 
  sll_classic_csl.F90 ../splines/contrib_rho.F90 ../splines/cubic_nonuniform_splines.F90)
target_link_libraries(vp_classic_csl
	sll_bgk_mesh
	sll_assertion 
        sll_constants 
        sll_file_io 
        sll_distribution_function 
        sll_mapped_meshes 
        sll_utilities 
        ode_solvers
        sll_scalar_field_2d 
        geometry_functions 
        sll_tridiagonal 
        sll_memory 
        sll_working_precision)

add_executable(vp_non_unif_mesh vp_non_unif_mesh.F90 
  sll_classic_csl.F90 ../splines/contrib_rho.F90 ../splines/cubic_nonuniform_splines.F90)
target_link_libraries(vp_non_unif_mesh
	sll_bgk_mesh
	sll_assertion 
        sll_constants 
        sll_file_io 
        sll_distribution_function 
        sll_mapped_meshes 
        sll_utilities 
        ode_solvers
        sll_scalar_field_2d 
        geometry_functions 
        sll_tridiagonal 
        sll_memory 
        sll_working_precision)


