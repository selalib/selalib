add_library( pppack STATIC
  banfac.f90 
  banslv.f90  
  bsplvb.f90
  woodbury.F90
  schur_complement.F90
)

install( TARGETS pppack ARCHIVE DESTINATION lib )

# Tests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif()
