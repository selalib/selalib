# Example (SLL_WARNING and SLL_ERROR macros)
ADD_EXECUTABLE( test_errors test_errors.F90 )
TARGET_LINK_LIBRARIES( test_errors sll_errors )

# [YG] Note: automatic test is missing because there is no way to tell CTest
#      that the unit test is successful if the executable exits abnormally.
