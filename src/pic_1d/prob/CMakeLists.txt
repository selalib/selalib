ADD_LIBRARY(probability_lgpl STATIC probability_lgpl.f90)

#PN Commented out to suppress a CMake warning
#TARGET_LINK_LIBRARIES( probability_lgpl  probability_lgpl )


#ADD_EXECUTABLE(test_probability_lgpl_prb  prob_prb.f90)
#ADD_DEPENDENCIES(   test_probability_lgpl_prb  probability_lgpl)
#TARGET_LINK_LIBRARIES( test_probability_lgpl_prb  probability_lgpl )
