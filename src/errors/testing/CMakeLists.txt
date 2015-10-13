# Example (SLL_WARNING and SLL_ERROR macros)
ADD_EXECUTABLE( test_crash test_crash.F90 )
TARGET_LINK_LIBRARIES( test_crash sll_errors )
