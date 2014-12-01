FIND_PROGRAM(LLRUN_EXECUTABLE
  NAMES llrun
  DOC "LoadLeveler tool to launch parallel jobs.")

IF (LLRUN_EXECUTABLE)
  SET (LLRUN_FOUND "YES")
ENDIF (LLRUN_EXECUTABLE)

MARK_AS_ADVANCED( LLRUN_FOUND LLRUN_EXECUTABLE)

IF(LLRUN_FOUND)

   SET(MPIEXEC_NUMPROC_FLAG "-procs") 
   MACRO(ADD_MPI_TEST TEST_NAME EXEC_NAME PROCS ARGS)
      ADD_TEST(NAME ${TEST_NAME}
               COMMAND "${LLRUN_EXECUTABLE} -k class=clallmds job_name=selalib total_tasks=${PROCS} node=1 node_usage=not_shared wall_clock_limit=00:10:00 job_type=mpich environment=COPY_ALL; mpirun " 

                       ${MPIEXEC_PREFLAGS}
	               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME}
	               ${MPIEXEC_POSTFLAGS} ${ARGS}
                       ${MPIEXEC_NUMPROC_FLAG} ${PROCS})
   ENDMACRO(ADD_MPI_TEST)

ELSE(LLRUN_FOUND)

   MACRO(ADD_MPI_TEST TEST_NAME EXEC_NAME PROCS ARGS)
      IF(NOT APPLE)
         STRING(REGEX REPLACE "mpiexec" "mpirun" MPIEXEC ${MPIEXEC})
      ENDIF()
      ADD_TEST(NAME ${TEST_NAME}
               COMMAND ${MPIEXEC} 
	               ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME}
	               ${MPIEXEC_POSTFLAGS} ${ARGS})
   ENDMACRO(ADD_MPI_TEST)

ENDIF(LLRUN_FOUND)


