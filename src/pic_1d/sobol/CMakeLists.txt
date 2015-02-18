ADD_LIBRARY(sobol_lgpl STATIC sobol.f90)



ADD_EXECUTABLE(test_sobol_prb  sobol_prb.f90)
ADD_DEPENDENCIES(   test_sobol_prb  sobol_lgpl)
TARGET_LINK_LIBRARIES( test_sobol_prb  sobol_lgpl )
