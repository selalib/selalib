add_library(sll_csl STATIC sll_csl.F90)
add_executable(test_csl unit_test.F90)
target_link_libraries(test_csl sll_csl sll_constants sll_splines 
  sll_distribution_function sll_mapped_meshes sll_utilities ode_solvers
  sll_scalar_field_2d geometry_functions sll_low_level_file_io sll_memory sll_assertion sll_working_precision)


add_library(sll_bgk_mesh STATIC sll_bgk_mesh.F90)

add_executable(vp_csl_bloc_unif vp_csl_bloc_unif.F90)
target_link_libraries(vp_csl_bloc_unif
	sll_bgk_mesh
	sll_csl
	sll_distribution_function
	sll_constants
	sll_splines
	sll_low_level_file_io 
	sll_mapped_meshes
	sll_utilities
	ode_solvers
	sll_scalar_field_2d
	geometry_functions
	sll_assertion sll_memory sll_working_precision)

#add_library(sll_classic_csl STATIC sll_classic_csl.F90)

add_executable(vp_classic_csl bgk_classic_csl.F90 
  sll_classic_csl.F90 ../splines/contrib_rho.F90 ../splines/cubic_nonuniform_splines.F90)
target_link_libraries(vp_classic_csl
	sll_bgk_mesh
	sll_assertion sll_constants sll_low_level_file_io 
  sll_distribution_function sll_mapped_meshes sll_utilities ode_solvers
  sll_scalar_field_2d geometry_functions sll_tridiagonal sll_memory sll_working_precision)

add_executable(vp_non_unif_mesh vp_non_unif_mesh.F90 
  sll_classic_csl.F90 ../splines/contrib_rho.F90 ../splines/cubic_nonuniform_splines.F90)
target_link_libraries(vp_non_unif_mesh
	sll_bgk_mesh
	sll_assertion sll_constants sll_low_level_file_io 
  sll_distribution_function sll_mapped_meshes sll_utilities ode_solvers
  sll_scalar_field_2d geometry_functions sll_tridiagonal sll_memory sll_working_precision)



#  use numeric_constants
#  use cubic_nonuniform_splines
#  use bgk_mesh_construction
#  use contrib_rho_module
#  use classical_conservative_semi_lagrangian



#env.Library('sll_memory',["../memory/sll_memory.F90"])
#env.Library('sll_precision', ["../precision/sll_working_precision.F90"])
#env.Library('sll_assert', ["../assert/sll_assert.F90"])
#env.Library('sll_ode_solvers',["../ode_solvers/sll_ode_solvers.F90"])
#env.Library('sll_tridiagonal', ["../numerical_utilities/tridiagonal.F90"])
#env.Library('sll_splines', ["../splines/sll_splines.F90","../splines/cubic_nonuniform_splines.F90" ])
#env.Library('sll_constants', ["../constants/numeric_constants.F90"])
#env.Library('sll_utilities',["../utilities/misc_utils.F90"])
#env.Library('sll_diagnostics',Glob("../diagnostics/sll_*.F90"))
#env.Library('sll_mesh_types', ["../mesh_types/sll_mesh_types.F90","../mesh_types/geometry_functions.F90"])
#env.Library('sll_distribution_function', ["../distribution_function/sll_distribution_function.F90"])



#env.Program('unit_test', ["unit_test.F90", "sll_csl.F90"])
#env.Program('unit_test_from_array', ["unit_test_from_array.F90", "sll_csl.F90"])
#env.Program('bgk_csl', ["bgk_csl.F90", "sll_csl.F90","sll_bgk_mesh.F90","../splines/contrib_rho.F90"])
#env.Program('bgk_classic_csl', ["bgk_classic_csl.F90", "sll_classic_csl.F90","sll_bgk_mesh.F90","../splines/contrib_rho.F90"])
#env.Program('bgk_csl_uniform', ["bgk_csl_uniform.F90", "sll_csl.F90","../splines/contrib_rho.F90"])
#env.Program('bgk_fv_unif', ["bgk_fv_unif.F90", "sll_csl.F90","../splines/contrib_rho.F90"])
#env.Program('bgk_fv', ["bgk_fv.F90", "sll_csl.F90","../splines/contrib_rho.F90"])
#env.Program('bgk_fv2', ["bgk_fv2.F90", "sll_csl.F90","../splines/contrib_rho.F90"])
#env.Program('bgk_cfd', ["bgk_CFD.F90", "sll_CFD.F90","../splines/contrib_rho.F90","sll_bgk_mesh.F90"])
#env.Program('cg_csl_uniform', ["cg_csl_uniform.F90", "cg_csl_uniform_module.F90"])
