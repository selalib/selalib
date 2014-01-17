#ADD_LIBRARY(mudpack STATIC cud2.f cud24.f cud24cr.f cud24sp.f cud2cr.f cud2sp.f cud3.f cud34.f cud34sp.f cud3cr.f cud3ln.f cud3pn.f cud3sp.f cudcom.f cuh2.f cuh24.f cuh24cr.f cuh2cr.f cuh3.f cuh34.f mud2.f mud24.f mud24cr.f mud24sp.f mud2cr.f mud2sa.f mud2sp.f mud3.f mud34.f mud34sp.f mud3cr.f mud3ln.f mud3pn.f mud3sa.f mud3sp.f mudcom.f muh2.f muh24.f muh24cr.f muh2cr.f muh3.f muh34.f resc2.f resc2cr.f resc2sp.f resc3.f resc3sp.f resm2.f resm2cr.f resm2sp.f resm3.f resm3cr.f resm3sp.f)


IF(MUDPACK_ENABLED)

   ADD_LIBRARY(mudpack STATIC mudcom.f mud2.f mud24.f 
                                    muh2.f muh24.f 
                                    mud24sp.f mud2sp.f
                                    mud24cr.f mud2cr.f 
                                    muh24cr.f muh2cr.f )

   MACRO(MUDPACK_TEST MUDPACK_TEST_NAME MUDPACK_TEST_COMMAND)
     ADD_EXECUTABLE( ${MUDPACK_TEST_COMMAND} "${MUDPACK_TEST_COMMAND}.f")
     TARGET_LINK_LIBRARIES(${MUDPACK_TEST_COMMAND} mudpack)
     ADD_TEST(NAME ${MUDPACK_TEST_NAME} 
           COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${MUDPACK_TEST_COMMAND} )
     SET(${MUDPACK_TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
   ENDMACRO(MUDPACK_TEST)

   MUDPACK_TEST(test_mud2    tmud2)
   MUDPACK_TEST(test_mud24   tmud24)
   MUDPACK_TEST(test_mud2sp  tmud2sp)
   MUDPACK_TEST(test_mud24sp tmud24sp)
   MUDPACK_TEST(test_mud2cr  tmud2cr)
   MUDPACK_TEST(test_mud24cr tmud24cr)
   MUDPACK_TEST(test_muh2    tmuh2)
   MUDPACK_TEST(test_muh24   tmuh24)
   MUDPACK_TEST(test_muh24cr tmuh24cr)
   MUDPACK_TEST(test_muh2cr  tmuh2cr)
   
   #ADD_EXECUTABLE( tmud2sa     tmud2sa.f)
   #ADD_EXECUTABLE( tmud3       tmud3.f)
   #ADD_EXECUTABLE( tmud34      tmud34.f)
   #ADD_EXECUTABLE( tmud34sp    tmud34sp.f)
   #ADD_EXECUTABLE( tmud3cr     tmud3cr.f)
   #ADD_EXECUTABLE( tmud3sa     tmud3sa.f)
   #ADD_EXECUTABLE( tmud3sp     tmud3sp.f)
   #ADD_EXECUTABLE( tmuh3       tmuh3.f)
   #ADD_EXECUTABLE( tmuh34      tmuh34.f)

ENDIF(MUDPACK_ENABLED)
