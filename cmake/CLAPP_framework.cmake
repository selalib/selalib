INCLUDE(ExternalProject)

# Convenience macro: find another library in the CLAPP framework
macro( find_clapp_library _NAME )
  find_package( ${_NAME} QUIET HINTS ${CLAPP_DIR}/cmake  $ENV{CLAPP_DIR}/cmake )
  if( ${_NAME}_FOUND )
    set( CLAPP_${_NAME} ON CACHE BOOL "Use pre-installed CLAPP/${_NAME} library" )
  endif()
  mark_as_advanced( ${_NAME}_DIR )
endmacro( find_clapp_library )

IF( CLAPP )

  SET( CLAPP_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "CLAPP installation directory")
  
  # Find all libraries in framework (assuming they have been installed)
  find_clapp_library( CLAPPIO )
  find_clapp_library( PLAF    )
  find_clapp_library( SPL     )
  find_clapp_library( DISCO   )
  find_clapp_library( FEMA    )
  find_clapp_library( GLT     )
  find_clapp_library( SPIGA   )
  find_clapp_library( SLLDJANGO )

  OPTION(BUILD_CLAPP "Build CLAPP framework from sources" OFF)

  IF(BUILD_CLAPP)

    EXTERNALPROJECT_ADD(
      clappio_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/clappio.git
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CLAPP_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)

    EXTERNALPROJECT_ADD(
      plaf_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/plaf.git
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CLAPP_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)

    EXTERNALPROJECT_ADD(
      spl_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/spl.git
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CLAPP_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)

    EXTERNALPROJECT_ADD(
      disco_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/disco.git
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CLAPP_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)

    EXTERNALPROJECT_ADD(
      glt_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/glt.git
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CLAPP_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)

    EXTERNALPROJECT_ADD(
      fema_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/fema.git
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CLAPP_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)

    ADD_DEPENDENCIES(plaf_lib  clappio_lib)
    ADD_DEPENDENCIES(spl_lib   plaf_lib)
    ADD_DEPENDENCIES(disco_lib spl_lib)
    ADD_DEPENDENCIES(glt_lib   disco_lib)
    ADD_DEPENDENCIES(fema_lib  glt_lib)

  ENDIF(BUILD_CLAPP)

ENDIF()
