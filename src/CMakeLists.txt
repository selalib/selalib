IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} 
       CACHE STRING "Choose the type of build, options are: Debug Release.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug 
       CACHE STRING "Choose the type of build, options are: Debug Release.")
ENDIF()

IF(OPTIONS_FILE)
   INCLUDE(${OPTIONS_FILE})
   MESSAGE("OPTIONS_FILE ARE READ")
ENDIF()

PROJECT(SeLaLib)

ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(OPENMP_ENABLED        OFF CACHE BOOL   "Enable use of OpenMP library")
SET(MPI_ENABLED           ON  CACHE BOOL   "Enable use of MPI library")
SET(HDF5_ENABLED          ON  CACHE BOOL   "Use HDF5 library for data output")
SET(HDF5_PARALLEL_ENABLED OFF CACHE BOOL   "Use Parallel HDF5")
SET(FFTW_ENABLED          ON  CACHE BOOL   "Enable use of FFTW library")
SET(PASTIX_ENABLED        OFF CACHE BOOL   "Build PasTiX solvers ")
SET(SUITESPARSE_ENABLED   OFF CACHE BOOL   "Enable UMFpack solvers")
SET(BUILD_SIMULATIONS     ON  CACHE BOOL   "Build selalib simulations")
SET(BUILD_PACKAGE         OFF CACHE BOOL   "Build selalib package")
SET(BUILD_PARALUTION      OFF CACHE BOOL   "Build paralution library")
SET(BUILD_JOREK           OFF CACHE BOOL   "Build jorek django solver")
SET(PROCS                 4   CACHE STRING "Number of MPI processes for tests")
SET(JOREK_GIT_BRANCH      selalib-branch CACHE STRING "Jorek branch to clone")

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  GET_FILENAME_COMPONENT(INSTALL_DIR "${CMAKE_SOURCE_DIR}/../usr" ABSOLUTE)
  SET(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH
    "Default installation directory" FORCE)
ENDIF()

FIND_PACKAGE(PythonInterp 3)
SET(PYTHON3_FOUND FALSE)
MESSAGE(STATUS "PYTHON_EXECUTABLE:${PYTHON_EXECUTABLE}")
MESSAGE(STATUS "PYTHON_VERSION_STRING:${PYTHON_VERSION_STRING}")
IF(PYTHON_EXECUTABLE AND PYTHON_VERSION_STRING MATCHES "3.*")
   STRING(COMPARE GREATER ${PYTHON_VERSION_STRING} "3.0.0" PYTHON3_FOUND)
ENDIF()
MESSAGE(STATUS "PYTHON3_FOUND:${PYTHON3_FOUND}")

# Specify the output for binary
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add new directory to find the FIND<PACKAGE>.cmake files. CMAKE_SOURCE_DIR is
# the directory where the project sources are located. CMAKE_MODULE_PATH is 
# the list of directories in which to search for CMake modules.
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules ${CMAKE_MODULE_PATH} )

IF(CMAKE_BUILD_TYPE MATCHES Debug)
   ADD_DEFINITIONS(-DDEBUG)
ENDIF()

# Define folder to be deleted by make clean
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/modules/)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/bin/)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_BINARY_DIR}/include/)

ADD_CUSTOM_TARGET( clear
   COMMAND find . -name *.mod -delete
   COMMAND find . -name *.a -delete
   COMMAND find . -name *.o -delete
   COMMAND find . -name *.cmake -delete
   COMMAND find . -name *.xmf -delete
   COMMAND find . -name *.h5 -delete
   COMMAND find . -name *.gnu -delete
   COMMAND find . -name *.dat -delete
   COMMAND rm -rf bin
   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
   COMMENT "Delete all .a, .mod, .o and folders /modules, /bin" VERBATIM)

INCLUDE(FortranConfig)
INCLUDE(LapackConfig)
INCLUDE(FFTConfig)
INCLUDE(HDF5Config)
INCLUDE(DoxygenConfig)
INCLUDE(CTest)

IF(MPI_ENABLED)
   INCLUDE(MPIConfig)
   INCLUDE(ProcessorCount)
   INCLUDE(POEConfig)
ENDIF()

IF(PASTIX_ENABLED)
   FIND_PACKAGE(SCOTCH)
   FIND_PACKAGE(PASTIX)
ELSE()
   IF(SUITESPARSE_ENABLED)
      FIND_PACKAGE(SUITESPARSE)
   ENDIF(SUITESPARSE_ENABLED)
ENDIF(PASTIX_ENABLED)

MACRO(SLL_TEST _NAME)
   ADD_EXECUTABLE( test_${_NAME} "unit_test_${_NAME}.F90")
   TARGET_LINK_LIBRARIES(test_${_NAME} ${ARGN})
   ADD_TEST(NAME ${_NAME} COMMAND test_${_NAME})
   SET(${_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
ENDMACRO(SLL_TEST)

EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE RESULT)
STRING(REGEX REPLACE "([0-9][0-9])/([0-9][0-9])/([0-9][0-9][0-9][0-9]).*" "\\3\\2\\1" 
TODAY ${RESULT})
STRING(REGEX REPLACE "([0-9][0-9][0-9][0-9])[0-9][0-9][0-9][0-9]" "\\1" 
CPACK_PACKAGE_VERSION_MAJOR ${TODAY})
STRING(REGEX REPLACE "[0-9][0-9][0-9][0-9]([0-9][0-9])[0-9][0-9]" "\\1" 
CPACK_PACKAGE_VERSION_MINOR ${TODAY} )
STRING(REGEX REPLACE "[0-9][0-9][0-9][0-9][0-9][0-9]([0-9][0-9])" "\\1" 
CPACK_PACKAGE_VERSION_PATCH ${TODAY} )

IF(BUILD_PACKAGE)
   SET(CMAKE_BUILD_TYPE Release
       CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
   INCLUDE(CPackConfig) 
ELSE()
   ENABLE_TESTING()
ENDIF(BUILD_PACKAGE)

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/modules)

ADD_SUBDIRECTORY(package)
ADD_SUBDIRECTORY(fftpack)
ADD_SUBDIRECTORY(mudpack)
ADD_SUBDIRECTORY(fishpack)
ADD_SUBDIRECTORY(memory)
ADD_SUBDIRECTORY(working_precision)
ADD_SUBDIRECTORY(assert)
ADD_SUBDIRECTORY(errors)
ADD_SUBDIRECTORY(constants)
ADD_SUBDIRECTORY(boundary_condition_descriptors)
ADD_SUBDIRECTORY(utilities)
ADD_SUBDIRECTORY(file_io)
ADD_SUBDIRECTORY(timer)
ADD_SUBDIRECTORY(splines)
ADD_SUBDIRECTORY(meshes)
ADD_SUBDIRECTORY(lagrange_interpolation)
ADD_SUBDIRECTORY(hermite_interpolation)
ADD_SUBDIRECTORY(integration)
ADD_SUBDIRECTORY(fft)
ADD_SUBDIRECTORY(periodic_interpolation)
ADD_SUBDIRECTORY(deboor_splines)
ADD_SUBDIRECTORY(reduction)

ADD_SUBDIRECTORY(pic_particle_types)
ADD_SUBDIRECTORY(random_deviate_generators)
ADD_SUBDIRECTORY(pic_particle_initializers)
ADD_SUBDIRECTORY(pic_accumulators)
ADD_SUBDIRECTORY(pic_utilities)
ADD_SUBDIRECTORY(field_accumulators)
ADD_SUBDIRECTORY(CAID_coordinate_transformations)
#ADD_SUBDIRECTORY(CAID_coordinate_transformations/circle_mp5)	
#ADD_SUBDIRECTORY(CAID_coordinate_transformations/mhdeq)
#ADD_SUBDIRECTORY(CAID_coordinate_transformations/identity)
ADD_SUBDIRECTORY(CAID_coordinate_transformations/square_4p_n10)

ADD_SUBDIRECTORY(interpolators)
ADD_SUBDIRECTORY(coordinate_transformations)
ADD_SUBDIRECTORY(characteristics)
ADD_SUBDIRECTORY(advection)
ADD_SUBDIRECTORY(fcisl)
ADD_SUBDIRECTORY(gyroaverage)
ADD_SUBDIRECTORY(maxwell_solvers)
ADD_SUBDIRECTORY(ampere_solvers)
ADD_SUBDIRECTORY(quasi_neutral_solvers)
ADD_SUBDIRECTORY(descriptors)

#Fortran version is checked in CMakeModules/FortranConfig.cmake
IF(F2008)
  ADD_SUBDIRECTORY( vector_space )
  ADD_SUBDIRECTORY( ode_integrators )
ENDIF(F2008)

# Container for all particle methods (mostly PIC)
ADD_SUBDIRECTORY( particle_methods )

# New output facilities for creating an XDMF database (WIP)
ADD_SUBDIRECTORY( xdmf )

IF(MPI_ENABLED)
   
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   ADD_SUBDIRECTORY(collective)
   ADD_SUBDIRECTORY(file_io_parallel)
   ADD_SUBDIRECTORY(remap)
   ADD_SUBDIRECTORY(point_to_point_communications)
   ADD_SUBDIRECTORY(distribution_function_multipatch)
   ADD_SUBDIRECTORY(poisson_solvers_parallel)
   ADD_SUBDIRECTORY(pif)
   #PN ADD_SUBDIRECTORY(maxwell_solvers_parallel)

ENDIF()

ADD_SUBDIRECTORY(parallel_utilities)
ADD_SUBDIRECTORY(fields)
ADD_SUBDIRECTORY(meshes_multipatch)
ADD_SUBDIRECTORY(coordinate_transformations_multipatch)
ADD_SUBDIRECTORY(fields_multipatch)
ADD_SUBDIRECTORY(AGMG)
ADD_SUBDIRECTORY(sparse_matrix_manager)
ADD_SUBDIRECTORY(general_coordinate_elliptic_solvers)
ADD_SUBDIRECTORY(general_coordinate_elliptic_solvers_multipatch)
ADD_SUBDIRECTORY(poisson_solvers)
ADD_SUBDIRECTORY(multigrid)
ADD_SUBDIRECTORY(visu_pic)
ADD_SUBDIRECTORY(time_solvers)
ADD_SUBDIRECTORY(distribution_function)
ADD_SUBDIRECTORY(ode_solvers)
ADD_SUBDIRECTORY(advection_field)
ADD_SUBDIRECTORY(mesh_calculus)

IF(BUILD_PACKAGE)
  MESSAGE(STATUS "SeLaLib will be installed in directory ${CMAKE_INSTALL_PREFIX}")
ENDIF(BUILD_PACKAGE)

IF(BUILD_SIMULATIONS)
  GET_FILENAME_COMPONENT( INSTALL_DIR_SIMS "${CMAKE_SOURCE_DIR}/../usr" ABSOLUTE )
  SET(CMAKE_INSTALL_PREFIX_SIMS ${INSTALL_DIR_SIMS} CACHE PATH 
    "simulations install dir")
  ADD_SUBDIRECTORY(simulations)
  ADD_SUBDIRECTORY(eulerian-finite-volume)
ENDIF(BUILD_SIMULATIONS)

IF(BUILD_PARALUTION)
  INCLUDE(ExternalProject)
  SET(PARALUTION_VERSION 1.0.0 CACHE STRING "Paralution version number")
  EXTERNALPROJECT_ADD( paralution
     URL  http://www.paralution.com/downloads/paralution-${PARALUTION_VERSION}.tar.gz
     CMAKE_ARGS -DSUPPORT_CUDA=OFF -DBUILD_EXAMPLES=OFF
     INSTALL_COMMAND ""
  )
ENDIF()


IF(BUILD_JOREK)

  SET(GIT_LOGIN $ENV{USER} CACHE STRING "INRIA GForge login")
  INCLUDE(ExternalProject)

  SET(JOREK_GIT_REPOSITORY git+ssh://${GIT_LOGIN}@scm.gforge.inria.fr/gitroot/jorek/jorek.git)

  IF(APPLE)
    SET(JOREK_CMAKE_ARGS -DBUILD_TESTING=OFF -DJOREK_DEBUG_TRACE=OFF -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DJOREK_SAVE_MATRIX=OFF)
  ELSE(APPLE)
    SET(JOREK_CMAKE_ARGS -DBUILD_TESTING=OFF -DJOREK_DEBUG_TRACE=OFF -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DJOREK_SAVE_MATRIX=OFF -DPETSC_DIR=${PETSC_DIR} -DPETSC_ARCH=${PETSC_ARCH})
  ENDIF(APPLE)

  IF(BUILD_PACKAGE)

    EXTERNALPROJECT_ADD(
      jorek
      GIT_REPOSITORY ${JOREK_GIT_REPOSITORY}
      GIT_TAG ${JOREK_GIT_BRANCH}
      CMAKE_ARGS ${JOREK_CMAKE_ARGS} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      )

  ELSE()

    EXTERNALPROJECT_ADD(
      jorek
      GIT_REPOSITORY ${JOREK_GIT_REPOSITORY}
      GIT_TAG ${JOREK_GIT_BRANCH}
      CMAKE_ARGS ${JOREK_CMAKE_ARGS}
      INSTALL_COMMAND "")

  ENDIF()
  
  FIND_PACKAGE(PETSc)
  ADD_SUBDIRECTORY(jorek)

ENDIF()

IF($ENV{PFUNIT})

  INCLUDE(ExternalProject)

  EXTERNALPROJECT_ADD(
    pFUnit
    PREFIX pFunit
    GIT_REPOSITORY git://git.code.sf.net/p/pfunit/code
    CMAKE_ARGS -DINSTALL_PATH=${CMAKE_BINARY_DIR})

ENDIF()

SET(BUILD_SPRNG OFF CACHE BOOL   "Build SPRNG libraries")
IF(BUILD_SPRNG)
  FIND_PACKAGE(MPI REQUIRED CXX)
  SET(SPRNG_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH  "CXX compiler for SPRNG")
  SET(SPRNG_LIBS "-lmpi_cxx" CACHE STRING  "LIBS configure argument for SPRNG")
  IF(APPLE)
    MESSAGE(STATUS "Warning: Build SPRNG library on mac is tricky")
  ENDIF(APPLE)
  INCLUDE(ExternalProject)
  SET(SPRNG_VERSION 5 CACHE STRING "SPRNG version number")
  EXTERNALPROJECT_ADD( sprng
     URL  http://www.sprng.org/Version5.0/sprng${SPRNG_VERSION}.tar.bz2
     SOURCE_DIR ${CMAKE_BINARY_DIR}/sprng
     BINARY_DIR ${CMAKE_BINARY_DIR}/sprng
     UPDATE_COMMAND autoreconf -if
     CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/sprng/configure --prefix=${CMAKE_BINARY_DIR}/sprng --with-mpi=yes CXX=${SPRNG_CXX_COMPILER} LIBS=${SPRNG_LIBS}
     BUILD_COMMAND ${MAKE} 
     INSTALL_COMMAND ""
  )
ENDIF()

MESSAGE(STATUS "############################################################")
MESSAGE(STATUS "####")
MESSAGE(STATUS "#### CMAKE_BUILD_TYPE:${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "#### Fortran_COMPILER_NAME:${Fortran_COMPILER_NAME}")
MESSAGE(STATUS "#### CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
MESSAGE(STATUS "#### Fortran FLAGS ")
IF(CMAKE_BUILD_TYPE MATCHES "Debug") 
   MESSAGE(STATUS "#### ${CMAKE_Fortran_FLAGS_DEBUG}")
ELSE()
   MESSAGE(STATUS "#### ${CMAKE_Fortran_FLAGS_RELEASE}")
ENDIF()
MESSAGE(STATUS "####")
MESSAGE(STATUS "############################################################")

ADD_CUSTOM_TARGET( clean_outputs
   COMMAND find . -name *.mtv -delete
   COMMAND find . -name *.gnu -delete
   COMMAND find . -name *.dat -delete
   COMMAND find . -name *.xmf -delete
   COMMAND find . -name *.h5  -delete
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   COMMENT "Delete all .gnu, .dat, .mtv, .h5, .xmf " VERBATIM)
