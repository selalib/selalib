FIND_PACKAGE(MPI REQUIRED CXX)
IF(APPLE)

  MESSAGE(STATUS "Warning: Build SPRNG library on mac is tricky")
  MESSAGE(STATUS "openmpi must be configured with cxx bindings")
  SET(SPRNG_CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/sprng/configure --prefix=${CMAKE_BINARY_DIR}/sprng --with-mpi=yes )

ELSE()

  SET(SPRNG_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH  "CXX compiler for SPRNG")
  SET(SPRNG_LIBS "-lmpi_cxx" CACHE STRING  "LIBS configure argument for SPRNG")
  SET(SPRNG_CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/sprng/configure --enable-silent-rules --prefix=${CMAKE_BINARY_DIR}/sprng --with-mpi=yes CXX=${SPRNG_CXX_COMPILER} LIBS=${SPRNG_LIBS} "CXXFLAGS=-w -g -O2")

ENDIF(APPLE)

INCLUDE(ExternalProject)
SET(SPRNG_VERSION 5 CACHE STRING "SPRNG version number")

EXTERNALPROJECT_ADD( sprng
   URL  http://www.sprng.org/Version5.0/sprng${SPRNG_VERSION}.tar.bz2
   SOURCE_DIR ${CMAKE_BINARY_DIR}/sprng
   BINARY_DIR ${CMAKE_BINARY_DIR}/sprng
   CONFIGURE_COMMAND ${SPRNG_CONFIGURE_COMMAND}
   BUILD_COMMAND ${MAKE}
   INSTALL_COMMAND ""
)
