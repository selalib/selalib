
ExternalProject_Get_Property(jorek SOURCE_DIR)
ExternalProject_Get_Property(jorek BINARY_DIR)

SET(CMAKE_Fortran_MODULE_DIRECTORY "${BINARY_DIR}/modules")

LINK_DIRECTORIES(${BINARY_DIR}/lib)

SET(JOREK_MODEL   poisson_2d)
ADD_EXECUTABLE(jorek_${JOREK_MODEL} 
   model_def.F90
   testcase.F90
   element_contribution.F90
   diagnostics_computation.F90
   model.F90
   sll_jorek.F90 
   unit_test_jorek.F90)
  
ADD_DEPENDENCIES(jorek_${JOREK_MODEL} jorek)
  
SET(JOREK_LIBRARIES  assembly 
                     basis 
                     blackbox 
                     boundary_conditions 
                     coordinates 
                     evaluator 
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
                     spm)
  
TARGET_LINK_LIBRARIES(jorek_${JOREK_MODEL} ${JOREK_LIBRARIES} ${PETSC_LIBRARIES})
  
SET(PROCS 1)
SET(PARAMETERS ${SOURCE_DIR}/models/${JOREK_MODEL}/parameter.txt)
SET(GEOMETRY  ${SOURCE_DIR}/geometries/square_periodic/n4x4)
ADD_TEST( NAME jorek
          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS}
                    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/jorek_${JOREK_MODEL} 
                       --parameters ${PARAMETERS}
                       --geometry   ${GEOMETRY} )

SET_TESTS_PROPERTIES( jorek PROPERTIES TIMEOUT 5 )
