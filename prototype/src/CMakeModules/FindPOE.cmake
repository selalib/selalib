FIND_PROGRAM(POE_EXECUTABLE
  NAMES poe
  DOC "IBM tool to launch parallel jobs.")

IF (POE_EXECUTABLE)
  SET (POE_FOUND "YES")
ENDIF (POE_EXECUTABLE)

MARK_AS_ADVANCED( POE_FOUND POE_EXECUTABLE)

IF(POE_FOUND)

   SET(HOSTLIST "${CMAKE_BINARY_DIR}/host.list")
   FILE(REMOVE ${HOSTLIST})
   SET(NPROCS 16)
   WHILE(NPROCS GREATER 0)
     FILE(APPEND ${HOSTLIST} "localhost\n")
     MATH(EXPR NPROCS "${NPROCS} - 1" ) 
   ENDWHILE( NPROCS GREATER 0)

   MACRO(SET_PROCS NPROCS)
      SET(MPIEXEC "poe")
      SET(MPIEXEC_NUMPROC_FLAG "")
      SET(PROCS "")
      SET(MPIEXEC_PREFLAGS "")
      SET(MPIEXEC_POSTFLAGS "-procs")
      SET(ARGS "${NPROCS}")
   ENDMACRO(SET_PROCS)

ELSE(POE_FOUND)

   MACRO(SET_PROCS NPROCS)
      SET(PROCS "${NPROCS}")
   ENDMACRO(SET_PROCS)

ENDIF(POE_FOUND)

