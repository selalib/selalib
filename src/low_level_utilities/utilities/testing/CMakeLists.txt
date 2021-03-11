add_executable(test_utilities test_utilities.F90)
target_link_libraries(test_utilities sll_utilities)
add_test(NAME utilities COMMAND test_utilities)

add_executable(test_fornberg test_fornberg.F90 module_test_fornberg.F90)
target_link_libraries(test_fornberg sll_utilities)
add_test(NAME fornberg COMMAND test_fornberg)

add_executable(test_int2string test_int2string.F90)
target_link_libraries(test_int2string sll_utilities)
