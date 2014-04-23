ADD_SUBDIRECTORY(sobol)
ADD_SUBDIRECTORY(prob)
ADD_LIBRARY(sll_pic_1d STATIC particle_loading.F90
			       pic_1d.F90
                              unit_test.F90                 )
                              
ADD_LIBRARY(sll_pic_1d_fem_solver STATIC bspline_fem_solver_1d.F90)


TARGET_LINK_LIBRARIES(sll_pic_1d_fem_solver  sll_constants
					    sll_working_precision
					sll_memory
					sll_logical_meshes
					sll_fft
					 sll_splines
					  sll_integration
					   sll_utilities
					   sll_collective)
					
TARGET_LINK_LIBRARIES(sll_pic_1d   sll_pic_1d_fem_solver
				      sll_constants
					sll_working_precision
					sll_memory
					sll_logical_meshes
					# sll_poisson_solvers
					sll_fft
					 sll_splines
					  sll_integration
					   sll_utilities
					   sll_interpolators
					   sll_collective
					   sobol_lgpl probability_lgpl sll_visu_pic )

ADD_EXECUTABLE(test_pic_1d  unit_test.F90)
ADD_DEPENDENCIES(   test_pic_1d  sll_pic_1d sobol_lgpl probability_lgpl sll_pic_1d_fem_solver)
TARGET_LINK_LIBRARIES(test_pic_1d  sll_pic_1d probability_lgpl  )


ADD_EXECUTABLE(test_bspline_fem_solver_1d  unit_test_bspline_fem_solver_1d.F90)
ADD_DEPENDENCIES(   test_bspline_fem_solver_1d  sll_utilities  sll_memory  sll_working_precision sll_pic_1d_fem_solver)
TARGET_LINK_LIBRARIES(test_bspline_fem_solver_1d  sll_pic_1d_fem_solver sll_file_io)

#ADD_EXECUTABLE(test_fft_bug  test_fft_bug.f90)
#TARGET_LINK_LIBRARIES(test_fft_bug  sll_fft sll_working_precision sll_memory sll_utilities)
