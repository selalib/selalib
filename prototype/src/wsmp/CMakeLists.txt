SET(WSMP_ENABLED OFF             CACHE BOOL "Use Watson Sparse Matrix Library ")
SET(WSMP_HOME    $ENV{WSMP_HOME} CACHE PATH "Watson Sparse Matrix location"    )

IF(WSMP_ENABLED)

   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   FIND_LIBRARY(WSMP_LIBRARY NAMES pwsmpEM64T-i13ibm
                             HINTS ${WSMP_HOME} 
                             PATH_SUFFIXES lib/mpi.ibm
                             DOC "Watson Sparse Matrix Library"
               )

   MESSAGE(STATUS "WSMP_LIBRARY:${WSMP_LIBRARY}")

   ADD_EXECUTABLE(pwssmp_ex1 pwssmp_ex1.f)
   TARGET_LINK_LIBRARIES(pwssmp_ex1 ${WSMP_LIBRARY} ${LAPACK_LIBRARIES})
   ADD_EXECUTABLE(pwssmp_ex2 pwssmp_ex2.f)
   TARGET_LINK_LIBRARIES(pwssmp_ex2 ${WSMP_LIBRARY} ${LAPACK_LIBRARIES})
   ADD_EXECUTABLE(pwssmp_ex3 pwssmp_ex3.f)
   TARGET_LINK_LIBRARIES(pwssmp_ex3 ${WSMP_LIBRARY} ${LAPACK_LIBRARIES})
ENDIF()


