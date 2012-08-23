<<<<<<< HEAD

=======
>>>>>>> origin/prototype-devel
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

IF(PASTIX_FOUND)
   include_directories(${PASTIX_INCLUDE_DIRS})
   add_executable(test_pastix fmurge.F90)
   add_executable(fsimple fsimple.F90 utils.F90)
   add_executable(fstep-by-step fstep-by-step.F90 utils.F90)

   target_link_libraries(test_pastix ${MURGE_LIBRARIES} 
                                     ${PASTIX_LIBRARIES} 
                                     ${BLAS_LIBRARIES}
                                     ${PTSCOTCH_LIBRARIES})

<<<<<<< HEAD
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
=======
   if(UNIX AND NOT APPLE)
       target_link_libraries(test_pastix rt)
>>>>>>> origin/prototype-devel
   endif()

ENDIF()
