add_executable( test_gaussian_deviate test_gaussian_deviate.F90 )
target_link_libraries( test_gaussian_deviate sll_random_deviate_generators)
add_test( NAME pic_gaussian_deviate COMMAND test_gaussian_deviate )
set_tests_properties( pic_gaussian_deviate PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

add_executable( test_hammersley test_hammersley.F90 )
target_link_libraries( test_hammersley sll_random_deviate_generators)
add_test( NAME pic_hammersley COMMAND test_hammersley )
set_tests_properties( pic_hammersley PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )


FIND_PACKAGE(GSL)

IF(FGSL_FOUND AND NOT PGI)
  ADD_EXECUTABLE( test_fgsl_rng mod_unit.f90 rng.f90)
  TARGET_INCLUDE_DIRECTORIES( test_fgsl_rng PRIVATE ${FGSL_INCLUDES})
  TARGET_LINK_LIBRARIES( test_fgsl_rng ${FGSL_LIBRARIES} )
ENDIF()
