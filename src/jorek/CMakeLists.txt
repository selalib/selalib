ExternalProject_Get_Property(jorek SOURCE_DIR)
ExternalProject_Get_Property(jorek BINARY_DIR)

SET(CMAKE_Fortran_MODULE_DIRECTORY ${BINARY_DIR}/modules)
INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
LINK_DIRECTORIES(${BINARY_DIR}/lib)

SET(JOREK_MODEL   poisson_2d)

ADD_LIBRARY(sll_jorek jorek_interface.F90 sll_jorek.F90)
ADD_DEPENDENCIES(sll_jorek jorek)

  
SET(JOREK_LIBRARIES  aggregate_field
                     assembly 
                     basis 
                     blackbox 
                     boundary_conditions 
                     dispatching
                     coordinates 
                     fem
                     evaluator 
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
  
TARGET_LINK_LIBRARIES(sll_jorek ${JOREK_LIBRARIES} 
                                sll_utilities
                                sll_file_io
                                ${PETSC_LIBRARIES})


ADD_EXECUTABLE(jorek_${JOREK_MODEL} unit_test_jorek.F90)
TARGET_LINK_LIBRARIES(jorek_${JOREK_MODEL} sll_jorek)
ADD_EXECUTABLE(jorek_landau_2d landau_2d.F90)
TARGET_LINK_LIBRARIES(jorek_landau_2d sll_interpolators sll_jorek)
  
SET(PROCS 4)
SET(PARAMETERS parameter.txt)
SET(GEOMETRY   ${SOURCE_DIR}/geometries/square/hbezier/n16x16)
SET(SOLVER     solver.txt)
IF(NOT APPLE)
  STRING(REGEX REPLACE "mpiexec" "mpirun" MPIEXEC ${MPIEXEC})
ENDIF()
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

