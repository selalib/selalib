ADD_LIBRARY(pic_vp_2d2v_cart_tau STATIC zone.F90 quietstart.F90 
  particules.F90 diagno.F90 particules_m6.F90 mpi_module.F90 )
TARGET_LINK_LIBRARIES(pic_vp_2d2v_cart_tau pic_vp_2d2v_cart_tau sll_poisson_solvers)

ADD_EXECUTABLE(sim_pic_vp_2d2v_cart_tau main.F90)
TARGET_LINK_LIBRARIES(sim_pic_vp_2d2v_cart_tau pic_vp_2d2v_cart_tau)

ADD_EXECUTABLE(test_interpolation_m4 test_interpolation_m4.F90)
TARGET_LINK_LIBRARIES(test_interpolation_m4 pic_vp_2d2v_cart_tau)

ADD_EXECUTABLE(test_high_order_interpolation test_high_order_interpolation.F90)
TARGET_LINK_LIBRARIES(test_high_order_interpolation pic_vp_2d2v_cart_tau
                                                    sll_particle_mesh_coupling
                                                    sll_particle_types)

#SET(ARGS "${CMAKE_CURRENT_SOURCE_DIR}/input/data.nml")
#ADD_TEST(NAME test_pic2d COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_pic2d ${ARGS})
#SET_TESTS_PROPERTIES(test_pic2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
