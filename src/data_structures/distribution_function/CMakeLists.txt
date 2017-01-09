SET( DISTRIBUTION_SOURCES_LIBRARY 
               sll_m_scalar_field_1d_old.F90
               sll_m_scalar_field_2d_base_old.F90
               sll_m_scalar_field_2d_old.F90
               sll_m_fdistribu4d_dk.F90
               sll_m_landau_2d_initializer.F90
               sll_m_tsi_2d_initializer.F90
               sll_m_distribution_function.F90
               sll_m_scalar_field_initializers_base.F90)

IF(PGI)

  FIND_PROGRAM(PYTHON_EXECUTABLE NAMES python)
  FIND_PROGRAM(CPP_EXECUTABLE NAMES cpp)
  SET(PREPROCESS_SCRIPT "${CMAKE_SOURCE_DIR}/python/cut_long_lines.py")
  
  SET(_FILE sll_m_initial_distribution)
  MESSAGE(STATUS "Cut long lines of file ${_FILE}.F90")

  EXECUTE_PROCESS(COMMAND ${CPP_EXECUTABLE} "-Iinclude" "-I${CMAKE_CURRENT_SOURCE_DIR}" "-E" "-w" "-P" 
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     INPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${_FILE}.F90"
                     OUTPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${_FILE}_pgi.F90")

  EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} ${PREPROCESS_SCRIPT} 
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  LIST(APPEND DISTRIBUTION_SOURCES_LIBRARY sll_m_initial_distribution_pgi.F90)

ELSE(PGI)

  LIST(APPEND DISTRIBUTION_SOURCES_LIBRARY sll_m_initial_distribution.F90)

ENDIF(PGI)

IF(MPI_ENABLED)

  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
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

