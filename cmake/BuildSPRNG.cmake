find_package(MPI REQUIRED CXX)

set(SPRNG_CONFIGURE_COMMAND
    ${CMAKE_BINARY_DIR}/sprng/configure --prefix=${CMAKE_BINARY_DIR}/sprng
    --with-fortran --with-mpi LIBS=-lmpi_cxx CXX=${CMAKE_CXX_COMPILER})

if(APPLE)

  message(STATUS "SPRNG : Warning Build the library on mac is tricky")
  message(STATUS "SPRNG : installed mpi must be configured with cxx bindings")
  message(
    STATUS
      "SPRNG : Use GNU c++ as compiler not clang, check with mpic++ -showme")
  set(SPRNG_PATCH_COMMAND patch -p0 <
                          ${CMAKE_CURRENT_SOURCE_DIR}/sprng/sprng-mac.patch)

else()

  message(STATUS "SPRNG : Build on fedora does not work")
  message(
    STATUS
      "SPRNG : You can try to configure and build it manually in ${CMAKE_BINARY_DIR}/sprng "
  )
  set(SPRNG_PATCH_COMMAND "")

endif(APPLE)

include(ExternalProject)
set(SPRNG_VERSION
    5
    CACHE STRING "SPRNG version number")

ExternalProject_Add(
  sprng
  URL http://www.sprng.org/Version5.0/sprng${SPRNG_VERSION}.tar.bz2
  SOURCE_DIR ${CMAKE_BINARY_DIR}/sprng
  BINARY_DIR ${CMAKE_BINARY_DIR}/sprng
  PATCH_COMMAND ${SPRNG_PATCH_COMMAND}
  CONFIGURE_COMMAND ${SPRNG_CONFIGURE_COMMAND}
  BUILD_COMMAND ${MAKE}
  BUILD_ALWAYS 0
  INSTALL_COMMAND "")
