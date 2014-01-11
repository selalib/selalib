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

SET(SLL_BUILD_PACKAGE OFF CACHE BOOL "Build library without test programs")
SET(MUDPACK_ENABLED   OFF CACHE BOOL "Build mudpack library")
SET(FISHPACK_ENABLED  OFF CACHE BOOL "Build fishpack library")
SET(FORTRANCL_ENABLED OFF CACHE BOOL "Enable OpenCL test programs")
SET(PASTIX_ENABLED    OFF CACHE BOOL "Build PasTiX solvers ")
SET(PETSC_ENABLED     OFF CACHE BOOL "Enable PETSc solvers")

IF(NOT SLL_BUILD_PACKAGE)
   SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../../usr" CACHE PATH "Installation directory of the library" FORCE)
ENDIF(NOT SLL_BUILD_PACKAGE)

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
   COMMAND rm -rf bin
   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
   COMMENT "Delete all .a, .mod, .o and folders /modules, /bin" VERBATIM)

SET(STDF95_ENABLED OFF CACHE BOOL "Use F95 norm")

INCLUDE(FortranConfig)
INCLUDE(LapackConfig)
INCLUDE(FFTConfig)
INCLUDE(MPIConfig)
INCLUDE(HDF5Config)

IF(MPI_MODULE_ENABLED)
   INCLUDE(ProcessorCount)
   INCLUDE(POEConfig)
ENDIF()


# include's .h files
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)

ADD_SUBDIRECTORY(fftpack)

IF(NOT SLL_BUILD_PACKAGE)

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
   ADD_SUBDIRECTORY(WENO)
   ADD_SUBDIRECTORY(integration)
   ADD_SUBDIRECTORY(fft)
   ADD_SUBDIRECTORY(periodic_interpolation)
   ADD_SUBDIRECTORY(deboor_splines)
   ADD_SUBDIRECTORY(reduction)
   IF(NOT STDF95)
     ADD_SUBDIRECTORY(coordinate_transformations)
     ADD_SUBDIRECTORY(characteristics)
     ADD_SUBDIRECTORY(advection)
   ENDIF(NOT STDF95)
   ADD_SUBDIRECTORY(interpolators)

   IF(MPI_MODULE_ENABLED)
   
      SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
      ADD_SUBDIRECTORY(collective)
      ADD_SUBDIRECTORY(remap)
      ADD_SUBDIRECTORY(parallel_utilities)
      ADD_SUBDIRECTORY(point_to_point_communications)


      IF(NOT STDF95)
         ADD_SUBDIRECTORY(parallel_array_initializers)
      ENDIF(NOT STDF95)

   ENDIF()
   
   ADD_SUBDIRECTORY(quasi_neutral_solvers)
   ADD_SUBDIRECTORY(electric_field_accumulators)
   
   ADD_SUBDIRECTORY(mudpack)
   ADD_SUBDIRECTORY(fishpack)
   ADD_SUBDIRECTORY(opencl)


   ADD_SUBDIRECTORY(sparse)
   ADD_SUBDIRECTORY(pastix)
   ADD_SUBDIRECTORY(wsmp)
   ADD_SUBDIRECTORY(petsc_solvers)
   ADD_SUBDIRECTORY(linear_solvers)

   IF(NOT STDF95)

      ADD_SUBDIRECTORY(fields)
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

   ENDIF(NOT STDF95)


   MESSAGE(STATUS "############################################################")
   MESSAGE(STATUS "####")
   MESSAGE(STATUS "#### CMAKE_BUILD_TYPE:${CMAKE_BUILD_TYPE}")
   MESSAGE(STATUS "#### CMAKE_Fortran_COMPILER:${CMAKE_Fortran_COMPILER}")
   MESSAGE(STATUS "#### Fortran FLAGS ")
   
   IF(CMAKE_BUILD_TYPE MATCHES "Debug") 
      MESSAGE(STATUS "#### ${CMAKE_Fortran_FLAGS_DEBUG}")
   ELSE()
      MESSAGE(STATUS "#### ${CMAKE_Fortran_FLAGS_RELEASE}")
   ENDIF()

   MESSAGE(STATUS "####")
   MESSAGE(STATUS "############################################################")

   IF(STDF95)

      MESSAGE(STATUS " ")
      MESSAGE(STATUS "####################  F95  ###################")
      MESSAGE(STATUS "## Directories not compiled in Fortran 95 : ##")
      MESSAGE(STATUS "## advection_field                          ##")
      MESSAGE(STATUS "## Linrood                                  ##")
      MESSAGE(STATUS "## BSL                                      ##")
      MESSAGE(STATUS "## coordinate_transformations               ##")
      MESSAGE(STATUS "## mapped_meshes                            ##")
      MESSAGE(STATUS "## parallel_array_initializers              ##")
      MESSAGE(STATUS "## CSL                                      ##")
      MESSAGE(STATUS "## distribution_function                    ##")
      MESSAGE(STATUS "## maxwel_solvers                           ##")
      MESSAGE(STATUS "## ode_solvers                              ##")
      MESSAGE(STATUS "## poisson_solvers                          ##")
      MESSAGE(STATUS "## time_solvers                             ##")
      MESSAGE(STATUS "## package                                  ##")
      MESSAGE(STATUS "## simulation                               ##")
      MESSAGE(STATUS "## visu_pic                                 ##")
      MESSAGE(STATUS "## fields (except geometry_functions.F90)   ##")
      MESSAGE(STATUS "##############################################")

   ENDIF(STDF95)

   INCLUDE(SequentialTests)

   IF(MPI_MODULE_ENABLED)

      INCLUDE(ParallelTests)

   ENDIF()


   INCLUDE(DoxygenConfig)

ENDIF(NOT SLL_BUILD_PACKAGE)

ADD_SUBDIRECTORY(package)

MESSAGE(STATUS "SeLaLib will be installed in directory ${CMAKE_INSTALL_PREFIX}")
