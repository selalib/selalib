ADD_LIBRARY(pic_vp_2d2v_cart_tau STATIC zone.F90 quietstart.F90 
  particules.F90 )
TARGET_LINK_LIBRARIES(pic_vp_2d2v_cart_tau sll_poisson_solvers sll_pic_visu)

ADD_EXECUTABLE(sim_pic_vp_2d2v_cart_tau main.F90)
TARGET_LINK_LIBRARIES(sim_pic_vp_2d2v_cart_tau pic_vp_2d2v_cart_tau)

ADD_EXECUTABLE(test_high_order_interpolation test_high_order_interpolation.F90)
TARGET_LINK_LIBRARIES(test_high_order_interpolation pic_vp_2d2v_cart_tau
                                                    sll_particle_mesh_coupling
                                                    sll_particle_groups)
