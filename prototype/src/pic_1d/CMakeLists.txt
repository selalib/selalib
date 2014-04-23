ADD_LIBRARY(sll_pic_1d STATIC particle_loading.F90
                              bspline_fem_solver_1d.F90
                              unit_test.F90                 )

TARGET_LINK_LIBRARIES(sll_pic_1d   sll_constants
					sll_working_precision
					sll_memory
					sll_logical_meshes
					# sll_poisson_solvers
					sll_fft
					 sll_splines
					  sll_integration
					   sll_utilities  )

ADD_EXECUTABLE(test_pic_1d  unit_test.F90)
TARGET_LINK_LIBRARIES(test_pic_1d  sll_pic_1d )
ADD_DEPENDENCIES(   test_pic_1d  sll_pic_1d)
