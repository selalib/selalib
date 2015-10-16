ADD_LIBRARY(sll_common_array_initializers STATIC 
            sll_m_common_array_initializers.F90 )

TARGET_LINK_LIBRARIES(sll_common_array_initializers
                      sll_assert
                      sll_working_precision
                      sll_constants)

IF(MPI_ENABLED)

  ADD_LIBRARY( sll_parallel_utilities STATIC sll_m_buffer_loader_utilities.F90 )
  
  TARGET_LINK_LIBRARIES( sll_parallel_utilities sll_remap sll_collective ) 
   

  ADD_LIBRARY(sll_parallel_array_utilities STATIC 
              sll_m_parallel_array_initializer.F90
              sll_m_common_array_initializers.F90 )

  TARGET_LINK_LIBRARIES(sll_parallel_array_utilities 
                        sll_remap
                        sll_coordinate_transformations
                        sll_file_io_parallel
                        sll_file_io
                        sll_utilities
                        sll_meshes)
  IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(testing)
  ENDIF(BUILD_TESTING)

ENDIF()
