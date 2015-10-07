ADD_LIBRARY( sll_pic_interface STATIC
  sll_pic_base.F90
  sll_remapped_pic_base.F90
  )

TARGET_LINK_LIBRARIES( sll_pic_interface sll_working_precision )

if( MPI_ENABLED )
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
  add_subdirectory( pic_1d_sheath )
  add_subdirectory( particle_types )
  add_subdirectory( particle_initializers )
  add_subdirectory( bsl_lt_pic_utilities )
endif( MPI_ENABLED )

