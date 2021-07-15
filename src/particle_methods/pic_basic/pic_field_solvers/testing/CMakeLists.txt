ADD_EXECUTABLE( test_pic_poisson_2d test_pic_poisson_2d.F90 )
TARGET_LINK_LIBRARIES( test_pic_poisson_2d sll_pic_field_solvers sll_particle_mesh_coupling sll_poisson_solvers )

ADD_TEST( NAME pic_poisson_2d COMMAND test_pic_poisson_2d )
SET_TESTS_PROPERTIES(pic_poisson_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
