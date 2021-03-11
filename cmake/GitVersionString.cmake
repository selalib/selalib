# Define the preprocessor macro SLL_GIT_VERSION, containing a string with the
# current git branch and git hash.  The macro can be used as code output for
# documentation purposes.
find_package(Git)
if(GIT_FOUND)
  exec_program(
    "git" ${CMAKE_CURRENT_SOURCE_DIR}
    ARGS "describe --all --long --dirty --tags"
    OUTPUT_VARIABLE SLL_GIT_RAW)
  string(REGEX MATCH "[^/]+$" SLL_GIT_VERSION ${SLL_GIT_RAW})
  message(STATUS "SLL_GIT_VERSION:${SLL_GIT_VERSION}")
  add_definitions(-DSLL_GIT_VERSION="${SLL_GIT_VERSION}")
endif()
