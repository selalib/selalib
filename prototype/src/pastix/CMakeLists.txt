IF(PASTIX_FOUND)

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
                                   ${SCOTCH_LIBRARIES}
                                   ${BLAS_LIBRARIES})

     MESSAGE(STATUS "BLAS_LIBRARIES:${BLAS_LIBRARIES}")
     TARGET_LINK_LIBRARIES(fsimple ${PASTIX_LIBRARIES} 
                                   ${PTSCOTCH_LIBRARIES}
                                   ${SCOTCH_LIBRARIES}
                                   ${BLAS_LIBRARIES})

     TARGET_LINK_LIBRARIES(pastix ${PASTIX_LIBRARIES} 
                                  ${PTSCOTCH_LIBRARIES}
                                  ${SCOTCH_LIBRARIES}
                                  ${BLAS_LIBRARIES})

     TARGET_LINK_LIBRARIES(fstep-by-step ${PASTIX_LIBRARIES} 
                                         ${PTSCOTCH_LIBRARIES}
                                         ${SCOTCH_LIBRARIES}
                                         ${BLAS_LIBRARIES})

     ADD_LIBRARY(sll_pastix sll_pastix.F90 sll_murge.F90 )

     TARGET_LINK_LIBRARIES(sll_pastix  sll_collective
                                       ${PASTIX_LIBRARIES} 
                                       ${MURGE_LIBRARIES} 
                                       ${PTSCOTCH_LIBRARIES}
                                       ${SCOTCH_LIBRARIES}
                                       ${BLAS_LIBRARIES})

     IF(UNIX AND NOT APPLE)

       TARGET_LINK_LIBRARIES(sll_pastix    rt)
       TARGET_LINK_LIBRARIES(fmurge        rt)
       TARGET_LINK_LIBRARIES(fsimple       rt)
       TARGET_LINK_LIBRARIES(fstep-by-step rt)

     ENDIF()

     ADD_EXECUTABLE(test_pastix unit_test.F90)
     IF(USE_MKL)
        TARGET_LINK_LIBRARIES(test_pastix sll_pastix) 
     ELSE()
        TARGET_LINK_LIBRARIES(test_pastix sll_pastix ${LAPACK_LIBRARIES}) 
     ENDIF(USE_MKL)


ENDIF()
