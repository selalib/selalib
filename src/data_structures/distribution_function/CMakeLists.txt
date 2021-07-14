set(DISTRIBUTION_SOURCES_LIBRARY
    sll_m_scalar_field_1d_old.F90
    sll_m_scalar_field_2d_base_old.F90
    sll_m_scalar_field_2d_old.F90
    sll_m_fdistribu4d_dk.F90
    sll_m_landau_2d_initializer.F90
    sll_m_tsi_2d_initializer.F90
    sll_m_distribution_function.F90
    sll_m_profile_functions.F90
    sll_m_scalar_field_initializers_base.F90)

list(APPEND DISTRIBUTION_SOURCES_LIBRARY sll_m_initial_distribution.F90)

if(MPI_ENABLED)

  add_library(
    sll_distribution_function STATIC
    ${DISTRIBUTION_SOURCES_LIBRARY}
    sll_m_distribution_function_initializer_4d.F90)

  target_link_libraries(
    sll_distribution_function sll_coordinate_transformations sll_meshes
    sll_file_io sll_collective sll_remap)
  if(BUILD_TESTING)
    add_subdirectory(testing)
  endif(BUILD_TESTING)

else()

  add_library(sll_distribution_function STATIC ${DISTRIBUTION_SOURCES_LIBRARY})

  target_link_libraries(sll_distribution_function
                        sll_coordinate_transformations sll_meshes sll_file_io)

endif(MPI_ENABLED)
