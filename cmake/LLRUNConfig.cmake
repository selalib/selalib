find_program(
  LLRUN_EXECUTABLE
  NAMES llrun
  DOC "LoadLeveler tool to launch parallel jobs.")

if(LLRUN_EXECUTABLE)
  set(LLRUN_FOUND "YES")
endif(LLRUN_EXECUTABLE)

mark_as_advanced(LLRUN_FOUND LLRUN_EXECUTABLE)

if(LLRUN_FOUND)

  set(MPIEXEC_NUMPROC_FLAG "-procs")
  macro(ADD_MPI_TEST TEST_NAME EXEC_NAME PROCS ARGS)
    add_test(
      NAME ${TEST_NAME}
      COMMAND
        "${LLRUN_EXECUTABLE} -k class=clallmds job_name=selalib total_tasks=${PROCS} node=1 node_usage=not_shared wall_clock_limit=00:10:00 job_type=mpich environment=COPY_ALL; mpirun "
        ${MPIEXEC_PREFLAGS} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME}
        ${MPIEXEC_POSTFLAGS} ${ARGS} ${MPIEXEC_NUMPROC_FLAG} ${PROCS})
  endmacro(ADD_MPI_TEST)

else(LLRUN_FOUND)

  macro(ADD_MPI_TEST TEST_NAME EXEC_NAME PROCS ARGS)
    if(NOT APPLE)
      string(REGEX REPLACE "mpiexec" "mpirun" MPIEXEC ${MPIEXEC})
    endif()
    add_test(
      NAME ${TEST_NAME}
      COMMAND
        ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME} ${MPIEXEC_POSTFLAGS}
        ${ARGS})
  endmacro(ADD_MPI_TEST)

endif(LLRUN_FOUND)
