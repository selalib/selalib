SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
add_executable(dk4d main_DK4D.F90
     utils_DK4D.F90
     input_DK4D.F90
     mesh_DK4D.F90
     magnetconf_DK4D.F90
     equilibrium_DK4D.F90
     fdistribu_DK4D.F90
     field3d_DK4D.F90
     elecpot_DK4D_hybrid.F90
     elecpot_DK4D_polar.F90
     QNsolver_DK4D_hybrid.F90
     QNsolver_DK4D_polar.F90
     vlasov_DK4D_hybrid.F90
     vlasov_DK4D_polar.F90
     diagnostics_DK4D.F90
     simulation_DK4D.F90)

target_link_libraries( dk4d  sll_simulation_base sll_quasi_neutral_solvers 
   sll_general_coordinate_elliptic_solvers sll_multipatch sll_advection
   sll_coordinate_transformations sll_fields sll_splines sll_remap 
   sll_collective sll_file_io)
   
   
