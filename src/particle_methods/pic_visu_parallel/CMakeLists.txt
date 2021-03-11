if(MPI_ENABLED)

  add_library(sll_pic_visu_parallel STATIC sll_m_pic_visu_parallel.F90)
  target_link_libraries(sll_pic_visu_parallel sll_utilities sll_file_io
                        sll_assert sll_pic_visu sll_collective)

  if(BUILD_TESTING)
    add_subdirectory(testing)
  endif(BUILD_TESTING)

endif(MPI_ENABLED)
