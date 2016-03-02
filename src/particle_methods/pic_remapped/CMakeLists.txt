ADD_LIBRARY( sll_remapped_pic STATIC
  sll_m_remapped_pic_base.F90
  )

TARGET_LINK_LIBRARIES( sll_remapped_pic
  sll_pic_opt2d_accumulators
  )

if( MPI_ENABLED )
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
  add_subdirectory( bsl_lt_pic )
  add_subdirectory( simple_pic )
  add_subdirectory( remapped_pic_utilities )
endif( MPI_ENABLED )

