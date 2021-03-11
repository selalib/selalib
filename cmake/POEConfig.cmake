find_program(
  POE_EXECUTABLE
  NAMES poe
  DOC "IBM tool to launch parallel jobs.")

if(POE_EXECUTABLE)
  set(POE_FOUND "YES")
endif(POE_EXECUTABLE)

mark_as_advanced(POE_FOUND POE_EXECUTABLE)

if(POE_FOUND)

  macro(HOSTLIST_CREATE)
    set(HOSTLIST "${CMAKE_CURRENT_BINARY_DIR}/host.list")
    set(NPROCS 16)
    while(NPROCS GREATER 0)
      file(APPEND ${HOSTLIST} "localhost\n")
      math(EXPR NPROCS "${NPROCS} - 1")
    endwhile(NPROCS GREATER 0)
  endmacro(HOSTLIST_CREATE)

  set(MPIEXEC ${POE_EXECUTABLE})
  set(MPIEXEC_NUMPROC_FLAG "-procs")
  macro(ADD_MPI_TEST TEST_NAME EXEC_NAME PROCS ARGS)
    hostlist_create()
    add_test(
      NAME ${TEST_NAME}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND
        ${MPIEXEC} ${MPIEXEC_PREFLAGS}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME} ${MPIEXEC_POSTFLAGS}
        ${ARGS} ${MPIEXEC_NUMPROC_FLAG} ${PROCS})
  endmacro(ADD_MPI_TEST)

else(POE_FOUND)

  macro(ADD_MPI_TEST TEST_NAME EXEC_NAME PROCS ARGS)
    set(MPIEXEC_PREFLAGS "--oversubscribe")
    string(REGEX REPLACE "mpiexec" "mpirun" MPIEXEC ${MPIEXEC})
    add_test(
      NAME ${TEST_NAME}
      COMMAND
        ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${PROCS}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME} ${MPIEXEC_POSTFLAGS}
        ${ARGS})
  endmacro(ADD_MPI_TEST)

endif(POE_FOUND)
