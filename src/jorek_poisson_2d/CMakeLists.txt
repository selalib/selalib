SET(GIT_LOGIN $ENV{USER} CACHE STRING "GForge login")

INCLUDE(ExternalProject)
EXTERNALPROJECT_ADD(
   jorek
   GIT_REPOSITORY git+ssh://${GIT_LOGIN}@scm.gforge.inria.fr/gitroot/jorek/jorek.git
   GIT_TAG global-devel
   CMAKE_ARGS -DBUILD_TESTING=OFF
   INSTALL_COMMAND ""
)

SET(JOREK_MODEL   poisson_2d)

################################
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
################################

FIND_PACKAGE(PETSc)

SET(CMAKE_Fortran_MODULE_DIRECTORY "${JOREK_BINARY_DIR}/modules")

LINK_DIRECTORIES(${JOREK_BINARY_DIR}/lib)

ExternalProject_Get_Property(jorek SOURCE_DIR)
ExternalProject_Get_Property(jorek BINARY_DIR)


ADD_EXECUTABLE(${JOREK_MODEL} 
   model_def.F90
   testcase.F90
   element_contribution.F90
   diagnostics_computation.F90
   model.F90
   solver.F90
   main.F90 )
  
ADD_DEPENDENCIES(${JOREK_MODEL} jorek)
  
SET(JOREK_LIBRARIES  assembly 
                     basis 
                     boundary_conditions 
                     coordinates 
                     evaluator 
                     graph 
                     jorek_context 
                     jorek_glob 
                     jorek_param 
                     matrix 
                     mesh
                     numbering 
                     output
                     spm)
  
TARGET_LINK_LIBRARIES(${JOREK_MODEL} ${JOREK_LIBRARIES} ${PETSC_LIBRARIES})
  
ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_BINARY_DIR}/parameter.txt
    COMMAND ${CMAKE_COMMAND} -E copy ${JOREK_SOURCE_DIR}/${JOREK_MODEL}/parameter.txt 
      ${CMAKE_BINARY_DIR}/parameter.txt
    MAIN_DEPENDENCY jorek
    COMMENT "copy \"parameter.txt\"")

