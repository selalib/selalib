SET( DISTRIBUTION_SOURCES_LIBRARY 
               sll_m_scalar_field_1d_old.F90
               sll_m_scalar_field_2d_base_old.F90
               sll_m_scalar_field_2d_old.F90
               sll_m_fdistribu4d_dk.F90
               sll_m_initial_distribution.F90
               sll_m_landau_2d_initializer.F90
               sll_m_tsi_2d_initializer.F90
               sll_m_distribution_function.F90
               sll_m_scalar_field_initializers_base.F90)

IF(MPI_ENABLED)

  ADD_LIBRARY(sll_distribution_function STATIC 
                 ${DISTRIBUTION_SOURCES_LIBRARY}
                 sll_m_distribution_function_initializer_4d.F90)
  
  TARGET_LINK_LIBRARIES(sll_distribution_function sll_coordinate_transformations
                                                  sll_meshes
                                                  sll_file_io
                                                  sll_collective
                                                  sll_remap )
  IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(testing)
  ENDIF(BUILD_TESTING)

ELSE()

  ADD_LIBRARY(sll_distribution_function STATIC ${DISTRIBUTION_SOURCES_LIBRARY})
  
  TARGET_LINK_LIBRARIES(sll_distribution_function sll_coordinate_transformations
                                                  sll_meshes
                                                  sll_file_io)

ENDIF(MPI_ENABLED)

