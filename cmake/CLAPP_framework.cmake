INCLUDE(ExternalProject)

# Convenience macro: find another library in the CLAPP framework
macro( find_clapp_library _NAME )
  find_package( ${_NAME} QUIET HINTS ${CLAPP_DIR}/cmake )
  if( ${_NAME}_FOUND )
    set( CLAPP_${_NAME} ON CACHE BOOL "Use pre-installed CLAPP/${_NAME} library" )
  endif()
  mark_as_advanced( ${_NAME}_DIR )
endmacro( find_clapp_library )

IF( CLAPP )

  SET( CLAPP_DIR ${CMAKE_BINARY_DIR} CACHE PATH "CLAPP installation directory")
  
  # Find all libraries in framework (assuming they have been installed)
  find_clapp_library( CLAPPIO )
  find_clapp_library( PLAF    )
  find_clapp_library( SPL     )
  find_clapp_library( DISCO   )
  find_clapp_library( FEMA    )
  find_clapp_library( SPIGA   )

  if( NOT CLAP_CLAPPIO)
    EXTERNALPROJECT_ADD(
      clappio_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/clappio.git
      UPDATE_DISCONNECTED 1
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)
  endif( NOT CLAP_CLAPPIO)
  
  if( NOT CLAP_PLAF)
    EXTERNALPROJECT_ADD(
      plaf_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/plaf.git
      UPDATE_DISCONNECTED 1
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)
    ADD_DEPENDENCIES(plaf_lib clappio_lib)
  endif( NOT CLAP_PLAF)

  
  if( NOT CLAP_SPL)
    EXTERNALPROJECT_ADD(
      spl_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/spl.git
      UPDATE_DISCONNECTED 1
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)
    ADD_DEPENDENCIES(spl_lib plaf_lib)
  endif( NOT CLAP_SPL)

  
  if( NOT CLAP_DISCO)
    EXTERNALPROJECT_ADD(
      disco_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/disco.git
      UPDATE_DISCONNECTED 1
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)
    
    ADD_DEPENDENCIES(disco_lib spl_lib)
  endif( NOT CLAP_DISCO)
  
  if( NOT CLAP_GLT)
    EXTERNALPROJECT_ADD(
      glt_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/glt.git
      UPDATE_DISCONNECTED 1
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)

    ADD_DEPENDENCIES(glt_lib disco_lib)
  endif( NOT CLAP_GLT)

  if( NOT CLAP_FEMA)
    EXTERNALPROJECT_ADD(
      fema_lib
      PREFIX clapp
      GIT_TAG devel
      GIT_REPOSITORY git@gitlab.mpcdf.mpg.de:clapp/fema.git
      UPDATE_DISCONNECTED 1
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_TESTING=OFF)

    ADD_DEPENDENCIES(fema_lib glt_lib)
  endif( NOT CLAP_FEMA)

ENDIF()
