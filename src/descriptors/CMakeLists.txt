ADD_LIBRARY(sll_descriptors STATIC 
            sll_descriptors.F90)

TARGET_LINK_LIBRARIES(sll_descriptors sll_working_precision
					sll_memory)



#   IF(BUILD_TESTING)
    ADD_EXECUTABLE( test_descriptors  unit_test.F90 )
    TARGET_LINK_LIBRARIES( test_descriptors sll_descriptors )
#   ENDIF(BUILD_TESTING)