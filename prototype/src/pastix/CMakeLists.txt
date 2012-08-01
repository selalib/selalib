set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

find_library(RT_LIBRARIES NAMES rt)

IF(PASTIX_FOUND)
   include_directories(${PASTIX_INCLUDE_DIRS})
   add_executable(test_pastix fmurge.F90)
   target_link_libraries(test_pastix ${MURGE_LIBRARIES} 
                                     ${PASTIX_LIBRARIES} 
                                     ${BLAS_LIBRARIES}
                                     ${PTSCOTCH_LIBRARIES}
                                     ${RT_LIBRARIES})

ENDIF()
