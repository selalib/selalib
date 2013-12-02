IF(PASTIX_ENABLED)

  IF( DEFINED ENV{SCOTCH_ROOT} )

     SET(SCOTCH_ROOT $ENV{SCOTCH_ROOT})
     SET(PTSCOTCH_ROOT $ENV{SCOTCH_ROOT})

  ENDIF()

  IF( DEFINED ENV{PASTIX_ROOT} )

     SET(PASTIX_ROOT $ENV{PASTIX_ROOT})
     SET(MURGE_ROOT $ENV{PASTIX_ROOT})

  ENDIF()

  FIND_PACKAGE(RT QUIET)
  FIND_PACKAGE(SCOTCH QUIET)
  FIND_PACKAGE(PTSCOTCH QUIET)
  FIND_PACKAGE(MURGE QUIET)
  FIND_PACKAGE(PASTIX QUIET)

  IF(PASTIX_FOUND AND PTSCOTCH_FOUND AND MURGE_FOUND AND SCOTCH_FOUND)

     MESSAGE(STATUS "PASTIX FOUND")
     INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})

     SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

     include_directories(${PASTIX_INCLUDE_DIRS})
     add_executable(test_pastix fmurge.F90)
     add_executable(fsimple fsimple.F90 utils.F90)
     add_executable(fstep-by-step fstep-by-step.F90 utils.F90)

     target_link_libraries(test_pastix ${MURGE_LIBRARIES} 
                                     ${PASTIX_LIBRARIES} 
                                     ${BLAS_LIBRARIES}
                                     ${PTSCOTCH_LIBRARIES})

     target_link_libraries(fsimple -lmatrix_driver ${PASTIX_LIBRARIES} 
                                     ${BLAS_LIBRARIES}
                                     ${PTSCOTCH_LIBRARIES})

     target_link_libraries(fstep-by-step -lmatrix_driver ${PASTIX_LIBRARIES} 
                                     ${BLAS_LIBRARIES}
                                     ${PTSCOTCH_LIBRARIES})

     if(UNIX AND NOT APPLE)
       target_link_libraries(test_pastix rt)
       target_link_libraries(fsimple rt)
       target_link_libraries(fstep-by-step rt)
     endif()

   ENDIF()

ENDIF(PASTIX_ENABLED)
