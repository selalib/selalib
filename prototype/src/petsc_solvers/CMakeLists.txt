IF(PETSC_ENABLED)

   FIND_PACKAGE(PETSc)

   IF(PETSC_FOUND)

      INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
      ADD_SUBDIRECTORY(petsc_solvers)

      IF(Fortran_COMPILER_NAME MATCHES "gfortran")
        SET(CMAKE_Fortran_FLAGS "-fopenmp")
        SET(OMP_LIBRARIES "-lgomp")
      ELSEIF(CMAKE_Fortran_COMPILER MATCHES "ifort")
        SET(CMAKE_Fortran_FLAGS "-openmp")
        SET(OMP_LIBRARIES "-openmp")
      ENDIF()

      ADD_LIBRARY(sll_petsc_solvers STATIC petsc_module.F90)
      ADD_EXECUTABLE(simple_test_petsc unit_test.F90)
      TARGET_LINK_LIBRARIES(simple_test_petsc sll_petsc_solvers ${PETSC_LIBRARIES} 
                                                          ${OMP_LIBRARIES})
      ADD_EXECUTABLE(cylinder cylinder.F90)
      TARGET_LINK_LIBRARIES(cylinder ${PETSC_LIBRARIES} )
      ADD_EXECUTABLE(cube cube.F90)
      TARGET_LINK_LIBRARIES(cube ${PETSC_LIBRARIES} )
      ADD_EXECUTABLE(laplacien2d_par laplacien2d_par.F90)
      TARGET_LINK_LIBRARIES(laplacien2d_par sll_petsc_solvers ${PETSC_LIBRARIES} )
      ADD_EXECUTABLE(laplacien2d_seq laplacien2d_seq.F90)
      TARGET_LINK_LIBRARIES(laplacien2d_seq sll_petsc_solvers ${PETSC_LIBRARIES} )
      ADD_EXECUTABLE(snes_example snes_example.F90)
      TARGET_LINK_LIBRARIES(snes_example sll_petsc_solvers ${PETSC_LIBRARIES} )
      ADD_EXECUTABLE(ts_matrix_free ts_matrix_free.F90)
      TARGET_LINK_LIBRARIES(ts_matrix_free sll_petsc_solvers ${PETSC_LIBRARIES} )
      SET(CMAKE_Fortran_FLAGS "")

   ENDIF(PETSC_FOUND)

ENDIF()
