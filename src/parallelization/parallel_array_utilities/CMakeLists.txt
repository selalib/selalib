ADD_LIBRARY(sll_parallel_array_utilities STATIC 
            sll_m_parallel_array_initializer.F90 )

TARGET_LINK_LIBRARIES(sll_parallel_array_utilities 
                      sll_remap
                      sll_coordinate_transformations
                      sll_file_io_parallel
                      sll_file_io
                      sll_utilities
                      sll_meshes
                      sll_common_array_initializers)


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)