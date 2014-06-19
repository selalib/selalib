IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} 
       CACHE STRING "Choose the type of build, options are: Debug Release.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug 
       CACHE STRING "Choose the type of build, options are: Debug Release.")
ENDIF()

IF(OPTIONS_FILE)
INCLUDE(${OPTIONS_FILE})
ENDIF()

PROJECT(SLL)

ENABLE_LANGUAGE(Fortran)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(MPI_MODULE_ENABLED    ON  CACHE BOOL "Enable use of MPI library")
SET(HDF5_ENABLED          ON  CACHE BOOL "Use HDF5 library for data output")
SET(HDF5_PARALLEL_ENABLED OFF CACHE BOOL "Use Parallel HDF5")
SET(FFTW_ENABLED          ON  CACHE BOOL "Enable use of FFTW library")
SET(BUILD_WITHOUT_TEST    OFF CACHE BOOL "Build library without test programs")
SET(MUDPACK_ENABLED       ON  CACHE BOOL "Build mudpack library")
SET(FISHPACK_ENABLED      OFF CACHE BOOL "Build fishpack library")
SET(FORTRANCL_ENABLED     OFF CACHE BOOL "Enable OpenCL test programs")
SET(PASTIX_ENABLED        OFF CACHE BOOL "Build PasTiX solvers ")
SET(PETSC_ENABLED         OFF CACHE BOOL "Enable PETSc solvers")
SET(SUITESPARSE_ENABLED   OFF CACHE BOOL "Enable UMFpack solvers")

IF(NOT BUILD_WITHOUT_TEST)
   SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../../usr" CACHE PATH "Installation directory of the library" FORCE)
ENDIF(NOT BUILD_WITHOUT_TEST)

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
INCLUDE(MPIConfig)
INCLUDE(HDF5Config)

IF(MPI_MODULE_ENABLED)
   INCLUDE(ProcessorCount)
   INCLUDE(POEConfig)
ENDIF()

IF(PASTIX_ENABLED)
  FIND_PACKAGE(RT QUIET)
  FIND_PACKAGE(SCOTCH QUIET)
  FIND_PACKAGE(PASTIX QUIET)
ENDIF(PASTIX_ENABLED)

IF(SUITESPARSE_ENABLED)
   IF(PASTIX_ENABLED)
   MESSAGE(SEND_ERROR "PASTIX IS ENABLED")
   ENDIF(PASTIX_ENABLED)
   FIND_PACKAGE(SUITESPARSE)
ENDIF(SUITESPARSE_ENABLED)

# include's .h files
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)

ADD_SUBDIRECTORY(package)
ADD_SUBDIRECTORY(fftpack)
ADD_SUBDIRECTORY(sparse_matrix_interface)

IF(NOT BUILD_WITHOUT_TEST)

   ENABLE_TESTING()

   ADD_SUBDIRECTORY(memory)
   ADD_SUBDIRECTORY(precision)
   ADD_SUBDIRECTORY(assert)
   ADD_SUBDIRECTORY(constants)
   ADD_SUBDIRECTORY(boundary_condition_descriptors)
   ADD_SUBDIRECTORY(utilities)
   ADD_SUBDIRECTORY(file_io)
   ADD_SUBDIRECTORY(timer)
   ADD_SUBDIRECTORY(non_linear_solvers)
   ADD_SUBDIRECTORY(splines)
   ADD_SUBDIRECTORY(logical_meshes)
   ADD_SUBDIRECTORY(lagrange_interpolation)
   ADD_SUBDIRECTORY(hermite_interpolation)
   ADD_SUBDIRECTORY(fcisl)
   ADD_SUBDIRECTORY(WENO)
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
   ADD_SUBDIRECTORY(gyroaverage)

   IF(MPI_MODULE_ENABLED)
   
      SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
      ADD_SUBDIRECTORY(collective)
      ADD_SUBDIRECTORY(remap)
      ADD_SUBDIRECTORY(parallel_utilities)
      ADD_SUBDIRECTORY(point_to_point_communications)
      ADD_SUBDIRECTORY(parallel_array_initializers)

   ENDIF()
   
   ADD_SUBDIRECTORY(quasi_neutral_solvers)
   ADD_SUBDIRECTORY(mudpack)
   ADD_SUBDIRECTORY(fishpack)
   ADD_SUBDIRECTORY(opencl)
   ADD_SUBDIRECTORY(sparse)
   ADD_SUBDIRECTORY(pastix)
   ADD_SUBDIRECTORY(wsmp)
   ADD_SUBDIRECTORY(petsc_solvers)
   ADD_SUBDIRECTORY(linear_solvers)
   ADD_SUBDIRECTORY(fields)
   ADD_SUBDIRECTORY(logical_meshes_multipatch)
   ADD_SUBDIRECTORY(coordinate_transformations_multipatch)
   ADD_SUBDIRECTORY(fields_multipatch)
   ADD_SUBDIRECTORY(discontinuous_galerkin) 
   ADD_SUBDIRECTORY(general_coordinate_elliptic_solvers)
   ADD_SUBDIRECTORY(poisson_solvers)
   ADD_SUBDIRECTORY(maxwell_solvers)
   ADD_SUBDIRECTORY(multigrid)
   ADD_SUBDIRECTORY(BSL)
   ADD_SUBDIRECTORY(CSL)
   ADD_SUBDIRECTORY(visu_pic)
   ADD_SUBDIRECTORY(mapped_meshes) # deprecated
   ADD_SUBDIRECTORY(time_solvers)
   ADD_SUBDIRECTORY(distribution_function)
   ADD_SUBDIRECTORY(ode_solvers)
   ADD_SUBDIRECTORY(advection_field)
   ADD_SUBDIRECTORY(LinRood)

   ADD_SUBDIRECTORY(mesh_calculus)
   ADD_SUBDIRECTORY(simulation)
   ADD_SUBDIRECTORY(eulerian-finite-volume)

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

   INCLUDE(SequentialTests)

   IF(MPI_MODULE_ENABLED)

      INCLUDE(ParallelTests)

   ENDIF()


   INCLUDE(DoxygenConfig)

ENDIF(NOT BUILD_WITHOUT_TEST)


MESSAGE(STATUS "SeLaLib will be installed in directory ${CMAKE_INSTALL_PREFIX}")
