SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
add_executable(vp4d main_VP4D.F90 utils_VP4D.F90 input_VP4D.F90 mesh_VP4D.F90 equilibrium_VP4D.F90 fdistribu_VP4D.F90 field2d_VP4D.F90 elecpot_Bsplines_VP4D.F90 QNsolver_Bsplines_VP4D.F90 vlasov_VP4D.F90 diagnostics_VP4D.F90 simulation_VP4D.F90)
target_link_libraries( vp4d ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS} )
