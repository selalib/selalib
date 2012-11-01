
FILE(COPY ${SELALIB_DIR}/check_compiler.F90 DESTINATION ${PROJECT_SOURCE_DIR} )
TRY_RUN( RUN_RESULT_VAR
	 COMPILE_RESULT_VAR
         ${CMAKE_BINARY_DIR}
         ${CMAKE_SOURCE_DIR}/check_compiler.F90)


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

IF(Fortran_COMPILER STREQUAL "GFORTRAN")
  add_definitions(-DGFORTRAN)
  set(CMAKE_Fortran_FLAGS "-g -Wall -cpp -pedantic -ffree-line-length-none -std=f2003 -fall-intrinsics -fbounds-check")
  add_definitions(-DDEBUG -DMPIF90)
ELSEIF(Fortran_COMPILER STREQUAL "INTEL")
	set(CMAKE_Fortran_FLAGS " ")
ELSEIF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qextname=flush -qthreaded -qhalt=e -qxlf2003=polymorphic")
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

SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
