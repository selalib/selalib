    ADD_EXECUTABLE( test_fcisl test_fcisl.F90 )
    TARGET_LINK_LIBRARIES( test_fcisl 
      sll_fcisl
      sll_meshes 
      sll_advection)

    ADD_EXECUTABLE( test_derivative_2d_oblic test_derivative_2d_oblic.F90 )
    TARGET_LINK_LIBRARIES( test_derivative_2d_oblic 
      sll_fcisl
      sll_meshes 
      sll_advection)


    ADD_EXECUTABLE( test_fcisl_toroidal test_fcisl_toroidal.F90 )
    TARGET_LINK_LIBRARIES( test_fcisl_toroidal 
      sll_fcisl
      sll_meshes 
      sll_advection)

    ADD_EXECUTABLE( test_aligned_derivative_2d aligned_derivative_2d.F90 )
    TARGET_LINK_LIBRARIES( test_aligned_derivative_2d 
      sll_fcisl
      sll_meshes 
      sll_advection)
