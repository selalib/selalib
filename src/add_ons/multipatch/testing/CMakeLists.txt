add_subdirectory(CAID_coordinate_transformations)
add_subdirectory(CAID_coordinate_transformations/square_4p_n10)

add_executable(test_distribution_function_4d_multipatch
               test_distribution_function_4d_multipatch.F90)

target_link_libraries(test_distribution_function_4d_multipatch sll_multipatch
                      sll_initial_profiles sll_timer sll_file_io)

if(PYTHON3_FOUND)
  set(ARGS " ")
  add_mpi_test(distribution_function_4d_multipatch
               test_distribution_function_4d_multipatch ${PROCS} ${ARGS})
  set_tests_properties(distribution_function_4d_multipatch
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endif(PYTHON3_FOUND)

set(INPUT_FILES_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/CAID_coordinate_transformations/square_4p_n10/)
file(GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml")
file(COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_meshes_multipatch test_meshes_multipatch.F90)
target_link_libraries(test_meshes_multipatch sll_multipatch)
add_test(NAME meshes_multipatch COMMAND test_meshes_multipatch)
set_tests_properties(meshes_multipatch PROPERTIES PASS_REGULAR_EXPRESSION
                                                  "PASSED")

add_executable(test_coordinate_transformation_multipatch_2d
               test_coordinate_transformation_multipatch_2d.F90)

target_link_libraries(test_coordinate_transformation_multipatch_2d
                      sll_multipatch sll_file_io)

if(PYTHON3_FOUND)
  add_test(NAME coordinate_transformation_multipatch_2d
           COMMAND test_coordinate_transformation_multipatch_2d)
  set_tests_properties(coordinate_transformation_multipatch_2d
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endif(PYTHON3_FOUND)

set(INPUT_FILES_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/CAID_coordinate_transformations/circle_mp5_pts12/
)
file(GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml")
file(COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_scalar_field_multipatch_2d
               test_scalar_field_multipatch_2d.F90)
target_link_libraries(test_scalar_field_multipatch_2d sll_multipatch)

if(PYTHON3_FOUND)
  add_test(NAME scalar_field_multipatch_2d
           COMMAND test_scalar_field_multipatch_2d)
  set_tests_properties(scalar_field_multipatch_2d
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endif(PYTHON3_FOUND)

set(INPUT_FILES_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/CAID_coordinate_transformations/square_4p_n10/)
file(GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml")
file(COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(
  test_general_coordinate_elliptic_solver_multipatch
  m_multipatch_helper_functions.F90
  test_general_coordinate_elliptic_solver_multipatch.F90)

target_link_libraries(test_general_coordinate_elliptic_solver_multipatch
                      sll_multipatch)

if(SUITESPARSE_ENABLED)

else()
  if(PYTHON3_FOUND)
    add_test(NAME general_coordinate_elliptic_solver_multipatch
             COMMAND test_general_coordinate_elliptic_solver_multipatch)
    set_tests_properties(general_coordinate_elliptic_solver_multipatch
                         PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
  endif(PYTHON3_FOUND)
endif(SUITESPARSE_ENABLED)

set(INPUT_FILES_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/CAID_coordinate_transformations/circle_mp5_pts12/
)
file(GLOB INPUT_FILES "${INPUT_FILES_DIR}/*.nml")
file(COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
