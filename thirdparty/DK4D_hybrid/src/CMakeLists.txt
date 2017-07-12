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
target_link_libraries( dk4d ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS} )
   
   
