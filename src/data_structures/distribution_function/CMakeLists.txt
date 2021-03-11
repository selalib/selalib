set(DISTRIBUTION_SOURCES_LIBRARY
    sll_m_scalar_field_1d_old.F90
    sll_m_scalar_field_2d_base_old.F90
    sll_m_scalar_field_2d_old.F90
    sll_m_fdistribu4d_dk.F90
    sll_m_landau_2d_initializer.F90
    sll_m_tsi_2d_initializer.F90
    sll_m_distribution_function.F90
    sll_m_scalar_field_initializers_base.F90)

if(PGI)

  find_program(PYTHON_EXECUTABLE NAMES python)
  find_program(FPP_EXECUTABLE NAMES gfortran)
  set(PREPROCESS_SCRIPT "${CMAKE_SOURCE_DIR}/python/cut_long_lines.py")

  set(_FILE sll_m_initial_distribution)
  message(STATUS "Cut long lines of file ${_FILE}.F90")

  execute_process(
    COMMAND ${FPP_EXECUTABLE} "-Iinclude" "-I${CMAKE_CURRENT_SOURCE_DIR}" "-E"
            "-w" "-P" "${CMAKE_CURRENT_SOURCE_DIR}/${_FILE}.F90"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${_FILE}_pgi.F90")

  execute_process(COMMAND ${PYTHON_EXECUTABLE} ${PREPROCESS_SCRIPT}
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  list(APPEND DISTRIBUTION_SOURCES_LIBRARY sll_m_initial_distribution_pgi.F90)

else(PGI)

  list(APPEND DISTRIBUTION_SOURCES_LIBRARY sll_m_initial_distribution.F90)

endif(PGI)

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
