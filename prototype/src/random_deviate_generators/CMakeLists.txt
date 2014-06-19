#

add_library(sll_random_deviate_generators STATIC gaussian.F90 hammersley.F90)

add_dependencies( sll_random_deviate_generators
                  sll_assertion 
		  sll_memory
		  sll_working_precision)

# add test
add_executable(test_random_deviate_generators unit_test.F90)
target_link_libraries( test_random_deviate_generators
                      sll_random_deviate_generators
                      sll_assertion
		      sll_memory
		      sll_working_precision )