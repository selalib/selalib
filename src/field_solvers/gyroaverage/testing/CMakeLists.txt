set(INPUT_FILES zeros_bessel.txt)
file(COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_gyroaverage_2d_polar_hermite
               test_gyroaverage_2d_polar_hermite.F90)

add_executable(test_gyroaverage_2d_polar_splines
               test_gyroaverage_2d_polar_splines.F90)

add_executable(test_gyroaverage_2d_polar_pade
               test_gyroaverage_2d_polar_pade.F90)

add_executable(test_qn_2d_polar_splines test_qn_2d_polar_splines.F90)

add_executable(test_qn_2d_polar_precompute test_qn_2d_polar_precompute.F90)

target_link_libraries(test_gyroaverage_2d_polar_hermite sll_gyroaverage)
target_link_libraries(test_gyroaverage_2d_polar_splines sll_gyroaverage)
target_link_libraries(test_gyroaverage_2d_polar_pade sll_gyroaverage)
target_link_libraries(test_qn_2d_polar_splines sll_gyroaverage)
target_link_libraries(test_qn_2d_polar_precompute sll_gyroaverage)

add_test(NAME gyroaverage_2d_polar_hermite
         COMMAND test_gyroaverage_2d_polar_hermite)
add_test(NAME gyroaverage_2d_polar_splines
         COMMAND test_gyroaverage_2d_polar_splines)
add_test(NAME gyroaverage_2d_polar_pade COMMAND test_gyroaverage_2d_polar_pade)
add_test(NAME qn_2d_polar_splines COMMAND test_qn_2d_polar_splines)
add_test(NAME qn_2d_polar_precompute COMMAND test_qn_2d_polar_precompute)

set_tests_properties(gyroaverage_2d_polar_hermite
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
set_tests_properties(gyroaverage_2d_polar_splines
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
set_tests_properties(gyroaverage_2d_polar_pade
                     PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
set_tests_properties(qn_2d_polar_splines PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")
set_tests_properties(qn_2d_polar_precompute PROPERTIES PASS_REGULAR_EXPRESSION
                                                       "PASSED")
