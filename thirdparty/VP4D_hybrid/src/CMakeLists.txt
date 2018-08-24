SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
INCLUDE_DIRECTORIES(${MPI_Fortran_MOD_DIR})
add_executable(vp4d main_VP4D.F90 utils_VP4D.F90 input_VP4D.F90 mesh_VP4D.F90 equilibrium_VP4D.F90 fdistribu_VP4D.F90 field2d_VP4D.F90 elecpot_Bsplines_VP4D.F90 QNsolver_Bsplines_VP4D.F90 vlasov_VP4D.F90 diagnostics_VP4D.F90 simulation_VP4D.F90)
target_link_libraries( vp4d  sll_simulation_base sll_quasi_neutral_solvers 
   sll_general_coordinate_elliptic_solvers sll_multipatch
   sll_coordinate_transformations sll_fields sll_splines sll_remap 
   sll_collective sll_file_io)
