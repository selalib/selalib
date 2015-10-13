ADD_LIBRARY( sll_errors STATIC sll_errors.F90 )

# Tests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )
