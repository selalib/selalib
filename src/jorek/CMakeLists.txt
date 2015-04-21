ExternalProject_Get_Property(jorek SOURCE_DIR)
ExternalProject_Get_Property(jorek BINARY_DIR)

SET(CMAKE_Fortran_MODULE_DIRECTORY ${BINARY_DIR}/modules)
INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
LINK_DIRECTORIES(${BINARY_DIR}/lib)

SET(JOREK_MODEL   poisson_2d)

ADD_EXECUTABLE(jorek_${JOREK_MODEL} sll_jorek.F90 
                                    unit_test_jorek.F90)
  
ADD_DEPENDENCIES(jorek_${JOREK_MODEL} jorek)
  
SET(JOREK_LIBRARIES  assembly 
                     basis 
                     blackbox 
                     boundary_conditions 
                     dispatching
                     coordinates 
                     evaluator 
                     fem
                     field
                     greenbox 
                     graph 
                     jorek_glob 
                     jorek_param 
                     linear_solver
                     matrix 
                     mesh
                     numbering 
                     output
                     quadratures
                     space
                     spm)
  
TARGET_LINK_LIBRARIES(jorek_${JOREK_MODEL} ${JOREK_LIBRARIES} ${PETSC_LIBRARIES})
  
SET(PROCS 4)
SET(PARAMETERS parameter.txt)
SET(GEOMETRY   ${SOURCE_DIR}/geometries/square_periodic/n16x16)
SET(SOLVER     solver.txt)
ADD_TEST( NAME jorek
          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/jorek_${JOREK_MODEL} 
              --parameters ${PARAMETERS}
              --solver     ${SOLVER} 
              --geometry   ${GEOMETRY} )

SET_TESTS_PROPERTIES( jorek PROPERTIES TIMEOUT 30 )

# ... Copy parameters file
FILE( COPY parameter.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE( COPY solver.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

