IF(PASTIX_ENABLED)

  FIND_PACKAGE(RT QUIET)
  FIND_PACKAGE(SCOTCH QUIET)
  FIND_PACKAGE(PASTIX QUIET)

  IF(PASTIX_FOUND AND PTSCOTCH_FOUND AND MURGE_FOUND AND SCOTCH_FOUND)

     INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})
     SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
     ADD_DEFINITIONS(-DDISTRIBUTED)

     INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})

     ADD_EXECUTABLE(fmurge fmurge.F90)
     ADD_EXECUTABLE(fsimple fsimple.F90 utils.F90)
     ADD_EXECUTABLE(pastix pastix.F90)
     ADD_EXECUTABLE(fstep-by-step fstep-by-step.F90 utils.F90)

     TARGET_LINK_LIBRARIES(fmurge  ${MURGE_LIBRARIES}
                                   ${PASTIX_LIBRARIES}
                                   ${PTSCOTCH_LIBRARIES}
                                   ${BLAS_LIBRARIES})

     TARGET_LINK_LIBRARIES(fsimple ${PASTIX_LIBRARIES} 
                                   ${PTSCOTCH_LIBRARIES}
                                   ${BLAS_LIBRARIES})

     TARGET_LINK_LIBRARIES(pastix ${PASTIX_LIBRARIES} 
                                  ${PTSCOTCH_LIBRARIES}
                                  ${BLAS_LIBRARIES})

     TARGET_LINK_LIBRARIES(fstep-by-step ${PASTIX_LIBRARIES} 
                                         ${PTSCOTCH_LIBRARIES}
                                         ${BLAS_LIBRARIES})

     ADD_LIBRARY(sll_pastix sll_pastix.F90 sll_murge.F90 )

     TARGET_LINK_LIBRARIES(sll_pastix  sll_collective
                                       ${MURGE_LIBRARIES} 
                                       ${PASTIX_LIBRARIES} 
                                       ${PTSCOTCH_LIBRARIES}
                                       ${BLAS_LIBRARIES})

     IF(UNIX AND NOT APPLE)

       TARGET_LINK_LIBRARIES(sll_pastix    rt)
       TARGET_LINK_LIBRARIES(fmurge        rt)
       TARGET_LINK_LIBRARIES(fsimple       rt)
       TARGET_LINK_LIBRARIES(fstep-by-step rt)

     ENDIF()

     ADD_EXECUTABLE(test_pastix unit_test.F90)
     TARGET_LINK_LIBRARIES(test_pastix sll_pastix)


   ENDIF()

ENDIF(PASTIX_ENABLED)
