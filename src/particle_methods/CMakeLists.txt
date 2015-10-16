ADD_LIBRARY( sll_pic_interface STATIC sll_m_pic_base.F90 )
TARGET_LINK_LIBRARIES( sll_pic_interface sll_working_precision )

if( MPI_ENABLED )
  set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
  add_subdirectory( pic_1d_sheath )
  add_subdirectory( pic_2d_standard )
  add_subdirectory( pic_basic )
  add_subdirectory( remapped_pic )
endif( MPI_ENABLED )

ADD_SUBDIRECTORY(visu_pic)
