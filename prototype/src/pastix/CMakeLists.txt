IF(PASTIX_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   include_directories(${PASTIX_INCLUDE_DIRS})
   add_executable(test_pastix fmurge.F90)
   target_link_libraries(test_pastix ${MURGE_LIBRARIES} 
                                     ${PASTIX_LIBRARIES} 
                                     ${PTSCOTCH_LIBRARIES}
                                     ${SCOTCH_LIBRARIES}
                                     ${BLAS_LIBRARIES})
ENDIF()
