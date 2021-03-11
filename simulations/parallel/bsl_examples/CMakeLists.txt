set(LIBS sll_interpolators sll_remap sll_file_io_parallel sll_file_io)

add_executable(advection_1d1d_cubic_splines advection_1d1d_cubic_splines.F90)
target_link_libraries(advection_1d1d_cubic_splines ${LIBS})

add_executable(advection_1d1d_quintic_splines
               advection_1d1d_quintic_splines.F90)
target_link_libraries(advection_1d1d_quintic_splines ${LIBS})

add_executable(advection_1d1d_arbitrary_degree_splines
               advection_1d1d_arbitrary_degree_splines.F90)
target_link_libraries(advection_1d1d_arbitrary_degree_splines ${LIBS})

add_executable(advection_2d_arbitrary_degree_splines
               advection_2d_arbitrary_degree_splines.F90)
target_link_libraries(advection_2d_arbitrary_degree_splines ${LIBS})

add_executable(polar_advection_arbitrary_degree_splines
               polar_advection_arbitrary_degree_splines.F90)
target_link_libraries(polar_advection_arbitrary_degree_splines
                      sll_coordinate_transformations ${LIBS})
