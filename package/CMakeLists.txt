# Create local variables
set( SLL_LIBS "" )  # Internal libraries (part of the SeLaLib project)
set( SLL_EXT  "" )  # External libraries (built with SeLaLib)
set( SLL_DEPS "" )  # Dependencies (pre-installed libraries)

# Extract names of all libraries built in selalib
get_property( library_targets GLOBAL PROPERTY LIBRARY_TARGETS )

# Distinguish between "internal" and "external" libraries
foreach( lib ${library_targets} )
  if( ${lib} MATCHES "sll_[a-zA-Z0-9]" )
    list( APPEND SLL_LIBS ${lib} )
  else()
    list( APPEND SLL_EXT ${lib} )
  endif()
endforeach()

# Extract system dependencies
foreach( lib ${library_targets} )
  get_target_property( deps ${lib} INTERFACE_LINK_LIBRARIES )
  if( deps )
    foreach( dep ${deps} )
      list( FIND library_targets ${dep} idx )
      if( ${idx} EQUAL -1 )
        list( APPEND SLL_DEPS ${dep} )
      endif()
    endforeach()
  endif()
endforeach()

# Remove duplicates without breaking dependency order
list( REVERSE SLL_DEPS )
list( REMOVE_DUPLICATES SLL_DEPS )
list( REVERSE SLL_DEPS )

#===============================================================================

# Print names of all sub-libraries in Selalib
message( STATUS "SELALIB LIBRARIES:" )
foreach( lib ${SLL_LIBS} )
  message( STATUS "  ${lib}" )
endforeach()

# Print names of built-in external 
message( STATUS "BUILT-IN EXTERNAL LIBRARIES:" )
foreach( lib ${SLL_EXT} )
  message( STATUS "  ${lib}" )
endforeach()

# Print all external dependencies
message( STATUS "EXTERNAL DEPENDENCIES:" )
foreach( lib ${SLL_DEPS} )
  message( STATUS "  ${lib}" )
endforeach()

#===============================================================================

# # Repackage SeLaLib into a single library (every time 'make' is invoked)
# # TODO: this can be handled entirely by CMake
# add_custom_target( selalib
#   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/make_package.py
#   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#   COMMENT "Repackage all static libraries into a single archive: 'libselalib.a'"
#   VERBATIM )
# add_dependencies( selalib ${SLL_LIBS} )

# CMake custom command: repackage SeLaLib into single archive "libselalib.a"
# NOTE: Command does not run if file is newer than its dependencies
# TODO: this can be handled entirely by CMake
set( SLL_LIB "${CMAKE_CURRENT_BINARY_DIR}/libselalib.a" )
add_custom_command( OUTPUT ${SLL_LIB}
  DEPENDS ${SLL_LIBS}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/make_package.py
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Repackage all static libraries into a single archive: 'libselalib.a'"
  VERBATIM )

# CMake target: update "libselalib.a" if obsolete, otherwise do nothing
# USAGE: make selalib-archive
add_custom_target( selalib-archive ALL DEPENDS ${SLL_LIB} )

# Wrap "libselalib.a" with CMake imported library "selalib", which depends on
# above "selalib-archive" target and therefore creates archive only if needed.
# USAGE: target_link_libraries( <MY_LIB> selalib ${SLL_EXT} ${SLL_DEPS} )
add_library( selalib STATIC IMPORTED GLOBAL )
set_target_properties( selalib PROPERTIES IMPORTED_LOCATION ${SLL_LIB} )
add_dependencies( selalib selalib-archive )

#===============================================================================

# Store cache variables (available to simulations)
set( SLL_LIB  ${SLL_LIB}  CACHE INTERNAL "Path to full SeLaLib library" )
set( SLL_EXT  ${SLL_EXT}  CACHE INTERNAL "List of built-in external libraries" )
set( SLL_DEPS ${SLL_DEPS} CACHE INTERNAL "List of external dependencies" )

# Find header files
file( GLOB_RECURSE SLL_MACROS ${PROJECT_SOURCE_DIR}/sll_*.h )

# Installation (library, headers, mod files)
install( FILES ${SLL_LIB}    DESTINATION lib )
install( FILES ${SLL_MACROS} DESTINATION "include/selalib" )
install( DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/"
  DESTINATION "include/selalib" FILES_MATCHING PATTERN "sll_m_*.mod" )
