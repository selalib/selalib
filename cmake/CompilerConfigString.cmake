# Macro to write out compiler and compiler flags in simulations
# author: Klaus Reuther
set(SLL_COMP "${Fortran_COMPILER_NAME} ${Fortran_COMPILER_VERSION}")
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
   set(SLL_COMP "${SLL_COMP} (${CMAKE_Fortran_FLAGS_DEBUG})")
ELSE()
   set(SLL_COMP "${SLL_COMP} (${CMAKE_Fortran_FLAGS_RELEASE})")
ENDIF()
MESSAGE( STATUS "SLL_COMP:${SLL_COMP}" )
add_definitions( -DSLL_COMP="${SLL_COMP}" )
