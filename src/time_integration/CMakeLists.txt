# Fortran version is checked in cmake/FortranConfig.cmake
if( F2008 )
  add_subdirectory( ode_integrators )
endif( F2008 )

add_subdirectory( splitting_methods )
add_subdirectory( characteristics )