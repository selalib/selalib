ADD_LIBRARY(sll_common_array_initializers STATIC 
            sll_common_array_initializers_module.F90 )

TARGET_LINK_LIBRARIES(sll_common_array_initializers
                      sll_assert
                      sll_working_precision
                      sll_constants)

IF(MPI_ENABLED)

  ADD_LIBRARY( sll_parallel_utilities STATIC buffer_loader_utilities.F90 )
  
  TARGET_LINK_LIBRARIES( sll_parallel_utilities sll_remap sll_collective ) 
   

  ADD_LIBRARY(sll_parallel_array_utilities STATIC 
              sll_parallel_array_initializer_module.F90
              sll_common_array_initializers_module.F90 )

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
