
# Determine how-to install the modules. CMAKE_BINARY_DIR is the directory
# in which the make command is invoked.
SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Add the modules directory to the list of include directories
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

GET_FILENAME_COMPONENT(Fortran_COMPILER_NAME "${CMAKE_Fortran_COMPILER}" NAME)

IF(CMAKE_SYSTEM_NAME MATCHES "BlueGene*")
##########################################################
# Try to determine the compiler
   SET(Fortran_COMPILER "IBM")
ELSE()
TRY_RUN( RUN_RESULT_VAR
	 COMPILE_RESULT_VAR
         ${CMAKE_BINARY_DIR}
         ${CMAKE_CURRENT_SOURCE_DIR}/check_compiler.F90
)
ENDIF()

# COMPILE_RESULT_VAR is set to true if try_run succeed
# RUN_RESULT_VAR is a string that represent the exit status
# message(STATUS "TRY_RUN_STATUS : ${COMPILE_RESULT_VAR}, EXIT_STATUS : ${RUN_RESULT_VAR}")


IF(COMPILE_RESULT_VAR)

   IF(${RUN_RESULT_VAR} STREQUAL 20 )

      MESSAGE(STATUS "COMPILER IS GNU, FORTRAN 2003 NOT SUPPORTED")
      SET(Fortran_COMPILER "GFORTRAN")
      SET(STDF95_ENABLED ON)
      ADD_DEFINITIONS(-DNOF03SUPPORT)

   ELSEIF(${RUN_RESULT_VAR} STREQUAL 21 )

      MESSAGE(STATUS "COMPILER IS GNU, FORTRAN 2003 SUPPORTED")
      SET(Fortran_COMPILER "GFORTRAN")

   ELSEIF(${RUN_RESULT_VAR} STREQUAL 30 )

      MESSAGE(STATUS "COMPILER IS INTEL, FORTRAN 2003 NOT SUPPORTED")
      SET(STDF95_ENABLED ON)
      SET(Fortran_COMPILER "INTEL")

   ELSEIF(${RUN_RESULT_VAR} STREQUAL 31 )

      MESSAGE(STATUS "COMPILER IS INTEL, FORTRAN 2003 SUPPORTED")
      SET(Fortran_COMPILER "INTEL")

   ENDIF()

ENDIF()

IF (${CMAKE_Fortran_COMPILER} MATCHES "ifort")
   EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "-v" OUTPUT_VARIABLE source_path)
   STRING(REGEX MATCH "1[0-9]\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
ELSE()
   EXEC_PROGRAM(${CMAKE_Fortran_COMPILER} ARGS "--version" OUTPUT_VARIABLE source_path)
   STRING(REGEX MATCH "4\\.[0-9]\\.[0-9]" Fortran_COMPILER_VERSION ${source_path})
ENDIF()
MESSAGE(STATUS "Fortran_COMPILER : ${Fortran_COMPILER}-${Fortran_COMPILER_VERSION}")


IF(Fortran_COMPILER STREQUAL "GFORTRAN")
   ADD_DEFINITIONS(-DGFORTRAN)
   SET(CMAKE_Fortran_FLAGS_RELEASE "-w -ffree-line-length-none -fall-intrinsics -O3 -fopenmp")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics -fbounds-check -fbacktrace -ffpe-trap=zero,overflow -O0")

ELSEIF(Fortran_COMPILER STREQUAL "INTEL")

   SET(CMAKE_Fortran_FLAGS_RELEASE "-nowarn -O3 -xHost -ip -openmp")
   SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all,noarg_temp_created -fpe0 -traceback -ftrapuv -fpic")

ELSEIF(Fortran_COMPILER MATCHES "IBM")

   SET(CMAKE_Fortran_FLAGS_DEBUG "-qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic")
   SET(CMAKE_Fortran_FLAGS_RELEASE "-qsmp=omp -qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic")

ELSE()

   MESSAGE(STATUS "NO KNOWN FORTRAN COMPILER FOUND")

ENDIF()


#STRING(FIND ${CMAKE_Fortran_FLAGS} "-std=f95" VAR)
#IF("${VAR}" STREQUAL "-1")
#  SET(STDF95_ENABLED OFF)
#ELSE()
IF(STDF95_ENABLED)
  SET(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -fmax-identifier-length=63")
  add_definitions(-DSTDF95)
  IF(${Fortran_COMPILER} STREQUAL "GFORTRAN")
     SET(CMAKE_Fortran_FLAGS_DEBUG "-pedantic -std=f95 -fmax-identifier-length=63 -g -Wall -cpp -ffree-line-length-none -fall-intrinsics -fbounds-check")
     SET(STDF95 YES)
     ADD_DEFINITIONS(-DSTDF95)
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(CLEAR CMAKE_Fortran_COMPILER)
