ADD_SUBDIRECTORY(sobol)
ADD_SUBDIRECTORY(prob)
ADD_LIBRARY(sll_pic_1d STATIC particle_loading.F90
			       pic_1d.F90
                              unit_test.F90  
                              pic_postprocessing.F90
                              particle_description.F90)
                              
ADD_LIBRARY(sll_pic_1d_fem_solver STATIC bspline_fem_solver_1d.F90)

TARGET_LINK_LIBRARIES(sll_pic_1d_fem_solver  sll_constants
					    sll_working_precision
					sll_memory
					sll_logical_meshes
					sll_boundary_condition_descriptors
					sll_fft
					 sll_splines
					  sll_integration
					   sll_utilities
					   sll_collective
					    sll_poisson_solvers)

					    
ADD_LIBRARY(sll_pic_1d_field_solver STATIC pic_1d_field_solver.F90
						         particle_description.F90)

TARGET_LINK_LIBRARIES(sll_pic_1d_field_solver  sll_constants
					    sll_working_precision
					sll_memory
					sll_logical_meshes
					sll_boundary_condition_descriptors
					sll_fft
					 sll_splines
					  sll_integration
					   sll_utilities
					   sll_collective
					    sll_poisson_solvers)
					    


					    
					    
					    
ADD_EXECUTABLE(test_bspline_fem_solver_1d  unit_test_bspline_fem_solver_1d.F90)
ADD_DEPENDENCIES(   test_bspline_fem_solver_1d  sll_utilities  sll_memory  sll_working_precision sll_pic_1d_fem_solver)
TARGET_LINK_LIBRARIES(test_bspline_fem_solver_1d  sll_pic_1d_fem_solver sll_file_io)
					   
					   
					   
					   

##IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)
					   
      TARGET_LINK_LIBRARIES(sll_pic_1d   sll_pic_1d_field_solver
				      sll_constants
					sll_working_precision
					sll_memory
					sll_logical_meshes
					sll_boundary_condition_descriptors
					# sll_poisson_solvers
					 sll_splines
					  sll_integration
					   sll_utilities
					   sll_interpolators
					   sll_collective
					 #  sll_file_io_parallel
					#   ${HDF5_LIBRARIES}
					   sobol_lgpl probability_lgpl sll_visu_pic sll_timer )


ADD_EXECUTABLE(test_pic_1d  unit_test.F90)
ADD_DEPENDENCIES(   test_pic_1d  sll_pic_1d sobol_lgpl probability_lgpl  sll_collective  )
TARGET_LINK_LIBRARIES(test_pic_1d  sll_pic_1d probability_lgpl  sll_collective )

## ENDIF()




#ADD_EXECUTABLE(test_fft_bug  test_fft_bug.f90)
#TARGET_LINK_LIBRARIES(test_fft_bug  sll_fft sll_working_precision sll_memory sll_utilities)
