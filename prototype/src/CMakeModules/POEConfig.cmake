IF(NOT DEFINED PROCESSOR_COUNT)
  # Unknown:
  SET(PROCESSOR_COUNT 0)

  # Linux:
  set(cpuinfo_file "/proc/cpuinfo")
  if(EXISTS "${cpuinfo_file}")
    file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
    list(LENGTH procs PROCESSOR_COUNT)
  endif()

  # Mac:
  if(APPLE)
    find_program(cmd_sys_pro "system_profiler")
    if(cmd_sys_pro)
      execute_process(COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info)
      string(REGEX REPLACE "^.*Total Number of Cores: ([0-9]+).*$" "\\1"
        PROCESSOR_COUNT "${info}")
    endif()
  endif()

  # Windows:
  if(WIN32)
    set(PROCESSOR_COUNT "$ENV{NUMBER_OF_PROCESSORS}")
  endif()

endif()

MESSAGE(STATUS "NUMBER OF PROCESSORS = ${PROCESSOR_COUNT}")
FIND_PROGRAM(POE_EXECUTABLE
  NAMES poe
  DOC "IBM tool to launch parallel jobs.")

IF (POE_EXECUTABLE)
  SET (POE_FOUND "YES")
ENDIF (POE_EXECUTABLE)

MARK_AS_ADVANCED( POE_FOUND POE_EXECUTABLE)

IF(POE_FOUND)

   SET(HOSTLIST "${CMAKE_CURRENT_BINARY_DIR}/host.list")
   SET(NPROCS 16)
   WHILE(NPROCS GREATER 0)
     FILE(APPEND ${HOSTLIST} "localhost\n")
     MATH(EXPR NPROCS "${NPROCS} - 1" ) 
   ENDWHILE( NPROCS GREATER 0)

   SET(MPIEXEC ${POE_EXECUTABLE}) 
   SET(MPIEXEC_NUMPROC_FLAG "-procs") 
   MACRO(ADD_MPI_TEST TEST_NAME EXEC_NAME PROCS ARGS)
      ADD_TEST(NAME ${TEST_NAME}
               COMMAND ${MPIEXEC} 
                       ${MPIEXEC_PREFLAGS}
	               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME}
	               ${MPIEXEC_POSTFLAGS} ${ARGS}
                       ${MPIEXEC_NUMPROC_FLAG} ${PROCS})
   ENDMACRO(ADD_MPI_TEST)

ELSE(POE_FOUND)

   MACRO(ADD_MPI_TEST TEST_NAME EXEC_NAME PROCS ARGS)
      ADD_TEST(NAME ${TEST_NAME}
               COMMAND ${MPIEXEC} 
	               ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME}
	               ${MPIEXEC_POSTFLAGS} ${ARGS})
   ENDMACRO(ADD_MPI_TEST)

ENDIF(POE_FOUND)
