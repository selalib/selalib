add_executable(test_compression test_compression.F90)
target_link_libraries(test_compression sll_compression)
add_test(NAME compression COMMAND test_compression)
# SET_TESTS_PROPERTIES( compression PROPERTIES TIMEOUT 20 )
