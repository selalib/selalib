
# Determine how-to install the modules. CMAKE_BINARY_DIR is the directory
# in which the make command is invoked.
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Add the modules directory to the list of include directories
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# DGFORTRAN is used to set the way in which numbers get converted to strings
# by the preprocessor. gfotran uses "x" while other preprocessors use the
# cpp #x. Defining the GFORTRAN flag chooses the first.
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
IF(Fortran_COMPILER_NAME STREQUAL "gfortran")
  message(STATUS "gfortran compiler")
  add_definitions(-DGFORTRAN)
  add_definitions(-DDEBUG)
ENDIF()
add_definitions(-DMPIF90)

##########################################################
# Try to determine the compiler
try_run( RUN_RESULT_VAR
	 COMPILE_RESULT_VAR
         ${CMAKE_BINARY_DIR}
         ${CMAKE_SOURCE_DIR}/check_compiler.F90
)

set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# COMPILE_RESULT_VAR is set to true if try_run succeed
# RUN_RESULT_VAR is a string that represent the exit status
# message(STATUS "TRY_RUN_STATUS : ${COMPILE_RESULT_VAR}, EXIT_STATUS : ${RUN_RESULT_VAR}")

set(NO_FORTRAN_2003 NO)
IF(COMPILE_RESULT_VAR)
	IF(${RUN_RESULT_VAR} STREQUAL 20 )
		message(STATUS "COMPILER IS GNU, FORTRAN 2003 NOT SUPPORTED")
		set(Fortran_COMPILER "GFORTRAN")
		set(NO_FORTRAN_2003 YES)
		add_definitions(-DNOF03SUPPORT)
	ELSEIF(${RUN_RESULT_VAR} STREQUAL 21 )
		message(STATUS "COMPILER IS GNU, FORTRAN 2003 SUPPORTED")
		set(Fortran_COMPILER "GFORTRAN")
	ELSEIF(${RUN_RESULT_VAR} STREQUAL 30 )
		message(STATUS "COMPILER IS INTEL, FORTRAN 2003 NOT SUPPORTED")
		set(NO_FORTRAN_2003 YES)
		set(Fortran_COMPILER "INTEL")
	ELSEIF(${RUN_RESULT_VAR} STREQUAL 31 )
		message(STATUS "COMPILER IS INTEL, FORTRAN 2003 SUPPORTED")
		set(Fortran_COMPILER "INTEL")
	ENDIF()
ELSE()
	message(STATUS "UNABLE TO DETERMINE WHICH COMPILER IS USED")
ENDIF()

##########################################################
# find out which compiler we are using.
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
IF(Fortran_COMPILER STREQUAL "GFORTRAN")
	set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics -fbounds-check")
ELSEIF(Fortran_COMPILER STREQUAL "INTEL")
	set(CMAKE_Fortran_FLAGS "-C")
        ADD_DEFINITIONS(-DINTEL)
ELSEIF(Fortran_COMPILER_NAME STREQUAL "xlf")
	set(CMAKE_Fortran_FLAGS "-qextname=flush -qthreaded -qhalt=e")
        set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ELSE()
	message(STATUS "NO KNOWN FORTRAN COMPILER FOUND")
ENDIF()

# this checks the existing collection of compiler flags to find if the
# fortran95 standard was requested. 
STRING(FIND ${CMAKE_Fortran_FLAGS} "-std=f95" VAR)
IF("${VAR}" STREQUAL "-1")
  set(STDF95 NO)
ELSE()
  set(STDF95 YES)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fmax-identifier-length=63")
  add_definitions(-DSTDF95)
ENDIF()

SET(STDF95_ENABLED OFF CACHE BOOL "Use F95 norm")
IF(${STDF95_ENABLED} AND ${Fortran_COMPILER} STREQUAL "GFORTRAN")
  SET(CMAKE_Fortran_FLAGS "-pedantic -std=f95 -fmax-identifier-length=63 -g -Wall -cpp -ffree-line-length-none -fall-intrinsics -fbounds-check")
  set(STDF95 YES)
  add_definitions(-DSTDF95)
ENDIF()

ADD_DEFINITIONS(-DMPIF90)
