 ADD_EXECUTABLE(test_curl_curl_solver test_curl_curl_solver.F90)
 TARGET_LINK_LIBRARIES(test_curl_curl_solver sll_curl_curl_solvers)

 ADD_TEST(NAME curl_curl_solver COMMAND test_curl_curl_solver)
 SET_TESTS_PROPERTIES(curl_curl_solver PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 500)

 
