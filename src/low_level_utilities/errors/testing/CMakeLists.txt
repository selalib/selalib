# Example (SLL_WARNING and SLL_ERROR macros)
add_executable(test_errors test_errors.F90)
target_link_libraries(test_errors sll_errors)

# [YG] Note: automatic test is missing because there is no way to tell CTest
# that the unit test is successful if the executable exits abnormally.
