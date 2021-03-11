if(PGI)

  find_program(PYTHON_EXECUTABLE NAMES python)
  find_program(FPP_EXECUTABLE NAMES gfortran)
  set(PREPROCESS_SCRIPT "${CMAKE_SOURCE_DIR}/python/cut_long_lines.py")

  set(_FILE sll_m_remapper)
  message(STATUS "Cut long lines of file ${_FILE}.F90")

  execute_process(
    COMMAND ${FPP_EXECUTABLE} "-Iinclude" "-I${CMAKE_CURRENT_SOURCE_DIR}" "-E"
            "-w" "-P" "${CMAKE_CURRENT_SOURCE_DIR}/${_FILE}.F90"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${_FILE}_pgi.F90")

  execute_process(COMMAND ${PYTHON_EXECUTABLE} ${PREPROCESS_SCRIPT}
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  add_library(sll_remap STATIC sll_m_remapper_pgi.F90)

else(PGI)

  add_library(sll_remap STATIC sll_m_remapper.F90)

endif(PGI)

target_link_libraries(sll_remap sll_collective sll_utilities)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
