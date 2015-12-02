# Fortran version is checked in cmake/FortranConfig.cmake
if( FULL_FORTRAN2003 )
  add_subdirectory( ode_integrators )
endif( FULL_FORTRAN2003 )

add_subdirectory( splitting_methods )
add_subdirectory( characteristics )
add_subdirectory( ode_solvers )
