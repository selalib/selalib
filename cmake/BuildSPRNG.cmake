FIND_PACKAGE(MPI REQUIRED CXX)
IF(APPLE)

  MESSAGE(STATUS "SPRNG : Warning Build the library on mac is tricky")
  MESSAGE(STATUS "SPRNG : openmpi must be configured with cxx bindings")
  SET(SPRNG_PATCH_COMMAND patch -p0 < ${CMAKE_CURRENT_SOURCE_DIR}/sprng/sprng-mac.patch)
  SET(SPRNG_CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/sprng/configure --prefix=${CMAKE_BINARY_DIR}/sprng --with-mpi  LIBS=-lmpi_cxx )

ELSE()

  SET(SPRNG_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH  "CXX compiler for SPRNG")
  SET(SPRNG_PATCH_COMMAND patch -p0 < ${CMAKE_CURRENT_SOURCE_DIR}/sprng/sprng-mac.patch)
  SET(SPRNG_CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/sprng/configure --prefix=${CMAKE_BINARY_DIR}/sprng --with-mpi )

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
   INSTALL_COMMAND ""
)
