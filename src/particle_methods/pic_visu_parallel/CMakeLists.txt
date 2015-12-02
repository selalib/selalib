IF (MPI_ENABLED)  
ADD_LIBRARY( sll_pic_visu_parallel STATIC sll_m_pic_visu_parallel.F90 )
TARGET_LINK_LIBRARIES( sll_pic_visu_parallel
  sll_utilities
  sll_file_io
  sll_assert
  sll_pic_visu
  sll_collective)
ENDIF (MPI_ENABLED) 
 
