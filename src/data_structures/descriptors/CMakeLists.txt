ADD_LIBRARY(sll_descriptors STATIC 
            sll_m_descriptors.F90)

TARGET_LINK_LIBRARIES(sll_descriptors sll_memory)


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
