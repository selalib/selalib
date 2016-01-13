# This option works with a separate 'make package' command
#ADD_CUSTOM_TARGET( library 
#  COMMAND ${CMAKE_BUILD_TOOL} -C ${CMAKE_BINARY_DIR} all
#  VERBATIM )
#
#ADD_CUSTOM_TARGET( package
#  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/make_package.py
#  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#  COMMENT "Repackage all static libraries into a single archive: 'libselalib.a'"
#  VERBATIM )
#
#ADD_DEPENDENCIES( package library )

# This option works every time 'make' is invoked
ADD_CUSTOM_TARGET( package ALL
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/make_package.py
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Repackage all static libraries into a single archive: 'libselalib.a'"
  VERBATIM )

# Installation
INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/libselalib.a DESTINATION lib )
