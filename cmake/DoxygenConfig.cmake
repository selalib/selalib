# for the documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)

  get_filename_component(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)
  set(DOXYGEN_OUTPUT_DIR ${CMAKE_BINARY_DIR}/doc)
  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
  message(STATUS "The documentation is in ${DOXYGEN_OUTPUT_DIR}")

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile-dev
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-dev @ONLY)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile-user
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-user @ONLY)

  add_custom_target(
    doc-dev
    COMMAND ${DOXYGEN_EXECUTABLE} -u ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-dev
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-dev
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen - developer version"
    VERBATIM)

  add_custom_target(
    doc-user
    COMMAND ${DOXYGEN_EXECUTABLE} -u ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-user
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-user
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen - user version"
    VERBATIM)

else(DOXYGEN_FOUND)
  message(STATUS "DOXYGEN NOT FOUND")
endif(DOXYGEN_FOUND)
