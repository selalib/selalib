# HACK: THIS SHOULD BE MADE BETTER AND AT A LOWER LEVEL. Problem: when hdf5 is
# present and is parallel, we need to use the MPI compiler wrapper and this
# should propagate to all its clients, even if these were envisioned as
# sequential modules. pic_visu is an example.

add_library(sll_pic_visu STATIC sll_m_pic_visu.F90 sll_m_pic_viewer.F90)

target_link_libraries(sll_pic_visu sll_coordinate_transformations sll_file_io
                      sll_utilities sll_file_io sll_assert)

if(BUILD_TESTING AND MPI_ENABLED)
  add_subdirectory(testing)
endif(BUILD_TESTING AND MPI_ENABLED)
