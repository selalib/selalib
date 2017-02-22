FIND_PACKAGE(MPI REQUIRED CXX)

SET(SPRNG_CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/sprng/configure --prefix=${CMAKE_BINARY_DIR}/sprng --with-fortran --with-mpi  LIBS=-lmpi_cxx CXX=${CMAKE_CXX_COMPILER} )

IF(APPLE)

  MESSAGE(STATUS "SPRNG : Warning Build the library on mac is tricky")
  MESSAGE(STATUS "SPRNG : installed mpi must be configured with cxx bindings")
  MESSAGE(STATUS "SPRNG : Use GNU c++ as compiler not clang, check with mpic++ -showme")
  SET(SPRNG_PATCH_COMMAND patch -p0 < ${CMAKE_CURRENT_SOURCE_DIR}/sprng/sprng-mac.patch)

ELSE()

  MESSAGE(STATUS "SPRNG : Build on fedora does not work")
  MESSAGE(STATUS "SPRNG : You can try to configure and build it manually in ${CMAKE_BINARY_DIR}/sprng ")
  SET(SPRNG_PATCH_COMMAND "")

ENDIF(APPLE)

INCLUDE(ExternalProject)
SET(SPRNG_VERSION 5 CACHE STRING "SPRNG version number")

EXTERNALPROJECT_ADD( sprng
   URL  http://www.sprng.org/Version5.0/sprng${SPRNG_VERSION}.tar.bz2
   SOURCE_DIR ${CMAKE_BINARY_DIR}/sprng
   BINARY_DIR ${CMAKE_BINARY_DIR}/sprng
   PATCH_COMMAND ${SPRNG_PATCH_COMMAND}
   CONFIGURE_COMMAND ${SPRNG_CONFIGURE_COMMAND}
   BUILD_COMMAND ${MAKE}
   BUILD_ALWAYS 0
   INSTALL_COMMAND ""
)
