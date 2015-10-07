SET(SLL_DIR ${CMAKE_SOURCE_DIR})

FILE(GLOB_RECURSE SLL_MACROS ${PROJECT_SOURCE_DIR}/*.h )

SET(SLL_SOURCE_FILES
    ${SLL_DIR}/memory/sll_memory.F90
    ${SLL_DIR}/working_precision/sll_working_precision.F90
    ${SLL_DIR}/assert/sll_assert.F90
    ${SLL_DIR}/advection/sll_advection_1d_periodic.F90
    ${SLL_DIR}/advection/sll_advection_1d_non_uniform_cubic_splines.F90
    ${SLL_DIR}/advection/sll_advection_1d_BSL.F90
    ${SLL_DIR}/advection/sll_advection_1d_CSL.F90
    ${SLL_DIR}/advection/sll_advection_1d_PSM.F90
    ${SLL_DIR}/advection/sll_advection_1d_base.F90
    ${SLL_DIR}/advection/sll_advection_2d_BSL.F90
    ${SLL_DIR}/advection/sll_advection_2d_tensor_product.F90
    ${SLL_DIR}/advection/sll_advection_2d_CSL.F90
    ${SLL_DIR}/advection/sll_advection_2d_base.F90
    ${SLL_DIR}/advection/sll_advection_1d_ampere.F90
    ${SLL_DIR}/advection/sll_advection_1d_spectral.F90
    ${SLL_DIR}/constants/sll_constants.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_explicit_euler.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_explicit_euler_conservative.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_trapezoid.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_trapezoid_conservative.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_base.F90
    ${SLL_DIR}/characteristics/sll_characteristics_2d_verlet.F90
    ${SLL_DIR}/characteristics/sll_characteristics_2d_explicit_euler.F90
    ${SLL_DIR}/characteristics/sll_characteristics_2d_explicit_euler_conservative.F90
    ${SLL_DIR}/characteristics/sll_characteristics_2d_base.F90
    ${SLL_DIR}/errors/sll_errors.F90
    ${SLL_DIR}/utilities/sll_utilities.F90
    ${SLL_DIR}/utilities/sll_primitives.F90
    ${SLL_DIR}/utilities/tridiagonal.F90
    ${SLL_DIR}/utilities/fornberg.F90
    ${SLL_DIR}/utilities/sll_toeplitz_penta_diagonal.F90
    ${SLL_DIR}/utilities/banfac.f90
    ${SLL_DIR}/utilities/banslv.f90
    ${SLL_DIR}/lagrange_interpolation/sll_lagrange_interpolation.F90
    ${SLL_DIR}/io/file_io/sll_binary_io.F90
    ${SLL_DIR}/io/file_io/sll_gnuplot.F90
    ${SLL_DIR}/io/file_io/sll_hdf5_io_serial.F90
    ${SLL_DIR}/io/file_io/sll_xdmf.F90
    ${SLL_DIR}/io/file_io/sll_ascii_io.F90
    ${SLL_DIR}/io/file_io/sll_xml_io.F90
    ${SLL_DIR}/io/file_io/sll_plotmtv.F90
    ${SLL_DIR}/io/file_io/sll_tri_mesh_xmf.F90
    ${SLL_DIR}/visu_pic/sll_visu_pic.F90
    ${SLL_DIR}/timer/sll_timer.F90
    ${SLL_DIR}/fields/sll_scalar_field_1d_base.F90
    ${SLL_DIR}/fields/sll_scalar_field_1d.F90
    ${SLL_DIR}/field_accumulators/sll_field_accumulators_2d.F90
    ${SLL_DIR}/fields/sll_scalar_field_2d_base.F90
    ${SLL_DIR}/fields/sll_scalar_field_2d.F90
    ${SLL_DIR}/field_solvers/general_coordinate_elliptic_solvers/sll_general_coordinate_elliptic_solver.F90
    ${SLL_DIR}/splines/sll_bsplines.F90
    ${SLL_DIR}/splines/sll_cubic_splines.F90
    ${SLL_DIR}/splines/sll_arbitrary_degree_splines.F90
    ${SLL_DIR}/splines/cubic_nonuniform_splines.F90
    ${SLL_DIR}/deboor_splines/sll_m_deboor_splines_1d.F90
    ${SLL_DIR}/deboor_splines/sll_m_deboor_splines_2d.F90
    ${SLL_DIR}/hermite_interpolation/sll_hermite_interpolation_2d.F90
    ${SLL_DIR}/interpolators/sll_interpolators_1d_base.F90
    ${SLL_DIR}/interpolators/sll_interpolators_2d_base.F90
    ${SLL_DIR}/interpolators/sll_cubic_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_cubic_spline_interpolator_1d_nonuniform.F90
    ${SLL_DIR}/interpolators/sll_cubic_spline_interpolator_2d.F90
    ${SLL_DIR}/interpolators/sll_arbitrary_degree_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_arbitrary_degree_spline_interpolator_2d.F90
    ${SLL_DIR}/interpolators/sll_periodic_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_lagrange_interpolator_1d.F90
    ${SLL_DIR}/field_solvers/maxwell_solvers/sll_maxwell_solvers_base.F90
    ${SLL_DIR}/field_solvers/maxwell_solvers/sll_maxwell_2d_fdtd.F90
    ${SLL_DIR}/periodic_interpolation/periodic_interp.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_1d_periodic.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_2d_polar.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_2d_periodic_fftpack.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_1d_base.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_1d_periodic_solver.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_1d_polar_solver.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_2d_polar_solver.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_2d_base.F90
    ${SLL_DIR}/interfaces/fft/sll_fft_utils.F90
    ${SLL_DIR}/integration/rectangle.F90
    ${SLL_DIR}/integration/trapz.F90
    ${SLL_DIR}/integration/gauss_legendre.F90
    ${SLL_DIR}/integration/gauss_lobatto.F90
    ${SLL_DIR}/simulations/simulations_sequential/sim_CG_polar/sll_polar_operators.F90
    ${SLL_DIR}/simulations/simulations_sequential/sim_CG_polar/sll_advection_polar.F90
    ${SLL_DIR}/interfaces/sparse_matrix_manager/partition.F90
    ${SLL_DIR}/interfaces/sparse_matrix_manager/sll_sparse_matrix_module.F90
    ${SLL_DIR}/time_solvers/sll_time_splitting_coeff.F90
    ${SLL_DIR}/boundary_condition_descriptors/sll_boundary_condition_descriptors.F90
    ${SLL_DIR}/meshes/sll_meshes_base.F90
    ${SLL_DIR}/meshes/sll_cartesian_meshes.F90
    ${SLL_DIR}/meshes/sll_triangular_meshes.F90
    ${SLL_DIR}/meshes/sll_hex_meshes.F90
    ${SLL_DIR}/fields/scalar_field_initializer_base.F90
    ${SLL_DIR}/distribution_function/sll_scalar_field_2d_old.F90
    ${SLL_DIR}/coordinate_transformations/sll_common_coordinate_transformations.F90
    ${SLL_DIR}/coordinate_transformations/sll_coordinate_transformations_2d_nurbs.F90
    ${SLL_DIR}/coordinate_transformations/sll_coordinate_transformation_2d_base.F90	
    ${SLL_DIR}/coordinate_transformations/sll_coordinate_transformations_2d.F90
    ${SLL_DIR}/simulations/multipatch/meshes_multipatch/sll_meshes_multipatch.F90
    ${SLL_DIR}/simulations/multipatch/coordinate_transformations_multipatch/sll_coordinate_transformation_multipatch.F90
    ${SLL_DIR}/simulations/multipatch/fields_multipatch/sll_scalar_field_2d_multipatch.F90
    ${SLL_DIR}/simulations/multipatch/general_coordinate_elliptic_solvers_multipatch/sll_general_coordinate_elliptic_solver_multipatch.F90
    ${SLL_DIR}/interfaces/sparse_matrix_manager/sll_sparse_matrix_mp_module.F90
    ${SLL_DIR}/mesh_calculus/sll_mesh_calculus_module.F90
    ${SLL_DIR}/distribution_function/distribution_function_initializer_4d.F90
    ${SLL_DIR}/distribution_function/sll_fdistribu4d_DK.F90
    ${SLL_DIR}/reduction/sll_reduction.F90
    ${SLL_DIR}/random_deviate_generators/gaussian.F90
    ${SLL_DIR}/random_deviate_generators/hammersley.F90
    ${SLL_DIR}/random_deviate_generators/moment_matching.F90
    ${SLL_DIR}/random_deviate_generators/sobol/sll_sobol.F90
    ${SLL_DIR}/random_deviate_generators/prob/sll_prob.F90
)

IF(F2008)
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/vector_space/sll_vector_space_base.F90
    ${SLL_DIR}/vector_space/sll_vector_space_real_arrays.F90
    ${SLL_DIR}/ode_integrators/sll_ode_integrator_base.F90
    ${SLL_DIR}/ode_integrators/rk_explicit.F90)
ENDIF(F2008)

IF(Fortran_COMPILER MATCHES "IBM")
ELSE()
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/utilities/sol.f
    ${SLL_DIR}/utilities/muls.f)
ENDIF()

IF(${FFT_LIB} MATCHES "FFTPACK")
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} ${SLL_DIR}/interfaces/fft/fftpack/sll_fft.F90)
ELSEIF(${FFT_LIB} MATCHES "FFTW")
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} ${SLL_DIR}/interfaces/fft/fftw3/sll_fft.F90) 
ELSE()
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} ${SLL_DIR}/interfaces/fft/sllfft/sll_fft.F90)
ENDIF()

IF(FFTW_ENABLED)
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/interfaces/fft/fftw3.F90
    ${SLL_DIR}/field_solvers/maxwell_solvers/sll_maxwell_2d_pstd.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_poisson_2d_periodic_fftw.F90
  )
ENDIF()

SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
  ${SLL_DIR}/field_solvers/poisson_solvers/sll_mudpack.F90
)

IF(SUITESPARSE_ENABLED)
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
    ${SLL_DIR}/interfaces/sparse_matrix_manager/spm_umfpack.F90
    ${SLL_DIR}/interfaces/sparse_matrix_manager/sll_sparse_matrix_umfpack_module.F90
  )
ENDIF(SUITESPARSE_ENABLED)

IF(MPI_ENABLED)
 SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
  ${SLL_DIR}/parallelization/collective/sll_collective.F90
  ${SLL_DIR}/parallelization/remap/sll_remap.F90
  ${SLL_DIR}/parallelization/point_to_point_communications/sll_point_to_point_comms.F90
  ${SLL_DIR}/field_solvers/poisson_solvers_parallel/sll_poisson_2d_periodic_cartesian_par.F90
  ${SLL_DIR}/field_solvers/poisson_solvers_parallel/sll_poisson_polar_parallel.F90
  ${SLL_DIR}/simulations/simulation_base_class.F90
  ${SLL_DIR}/parallelization/parallel_utilities/sll_common_array_initializers_module.F90
  ${SLL_DIR}/parallelization/parallel_utilities/sll_parallel_array_initializer_module.F90
  ${SLL_DIR}/parallelization/parallel_utilities/buffer_loader_utilities.F90
  ${SLL_DIR}/io/file_io_parallel/sll_gnuplot_parallel.F90)
ENDIF(MPI_ENABLED)

IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

 SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
  ${SLL_DIR}/io/file_io_parallel/sll_hdf5_io_parallel.F90
  ${SLL_DIR}/io/file_io_parallel/sll_xdmf_parallel.F90
  ${SLL_DIR}/simulations/multipatch/distribution_function_multipatch/sll_distribution_function_4d_multipatch.F90
  ${SLL_DIR}/simulations/simulations_parallel/sim_4d_qns_general/simulation_4d_qns_general.F90
  ${SLL_DIR}/simulations/simulations_parallel/sim_2d_vlasov_poisson_cartesian/simulation_2d_vlasov_poisson_cartesian.F90
  ${SLL_DIR}/simulations/simulations_parallel/sl_1d1v_va_cart_spectral/simulation_2d_vlasov_ampere_cartesian.F90
  ${SLL_DIR}/simulations/multipatch/sim_4d_qns_general_multipatch/simulation_4d_qns_general_multipatch.F90
  ${SLL_DIR}/simulations/simulations_parallel/sim_4d_drift_kinetic_hybrid/simulation_4d_DK_hybrid.F90
  ${SLL_DIR}/simulations/simulations_parallel/sim_4d_drift_kinetic_polar/simulation_4d_drift_kinetic_polar.F90)

ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

IF(BUILD_PACKAGE)

 SET(BUILD_SHARED_LIBRARY  OFF CACHE BOOL "build shared library")
  IF(BUILD_SHARED_LIBRARY)
     ADD_LIBRARY(selalib SHARED ${SLL_SOURCE_FILES})
     TARGET_LINK_LIBRARIES(selalib dfftpack mudpack 
       ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${FFTW_LIBRARIES})
  ELSE()
     ADD_LIBRARY(selalib STATIC ${SLL_SOURCE_FILES})
     INSTALL(TARGETS selalib ARCHIVE DESTINATION lib COMPONENT libraries)
  ENDIF()
  
  
  SET(DIST_NAME "selalib-${CPACK_PACKAGE_VERSION}")
  
  ADD_CUSTOM_TARGET(release 
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
          ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/${DIST_NAME}/include
  COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_CURRENT_SOURCE_DIR}/macports/README 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/
  COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_CURRENT_SOURCE_DIR}/macports/CMakeLists.txt 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/
  COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_CURRENT_SOURCE_DIR}/macports/Doxyfile.in 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/Doxyfile
  COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_CURRENT_SOURCE_DIR}/../src/documentation.f90 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/
  COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_CURRENT_SOURCE_DIR}/selalib-mpi.h
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/include
  COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_CURRENT_SOURCE_DIR}/selalib.h
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/include
  COMMAND ${CMAKE_COMMAND} -E copy 
          ${CMAKE_CURRENT_SOURCE_DIR}/../doc/sphinx/_static/logo2.png 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
          ${CMAKE_CURRENT_SOURCE_DIR}/../src/external/fftpack 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/fftpack
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
          ${CMAKE_CURRENT_SOURCE_DIR}/../src/external/mudpack 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/mudpack
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
          ${CMAKE_CURRENT_SOURCE_DIR}/examples 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/examples
  COMMAND ${CMAKE_COMMAND} -E remove_directory 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/CMakeCache.txt
  COMMAND ${CMAKE_COMMAND} -E tar zcf ${DIST_NAME}.tar.gz 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}
  COMMAND ${CMAKE_COMMAND} -E md5sum ${DIST_NAME}.tar.gz
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/${DIST_NAME}
           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
           COMMENT "Generating Release" VERBATIM)
   
  ADD_CUSTOM_TARGET(make_source_dir 
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/${DIST_NAME}/src)
  ADD_DEPENDENCIES(release make_source_dir)
  
  FOREACH(_file ${SLL_SOURCE_FILES})
     get_filename_component(_target ${_file} NAME_WE)
     ADD_CUSTOM_TARGET(${_target}_f90 ${CMAKE_COMMAND} 
     -E copy ${_file} ${CMAKE_BINARY_DIR}/${DIST_NAME}/src)
     ADD_DEPENDENCIES(release ${_target}_f90)
     ADD_DEPENDENCIES(${_target}_f90 make_source_dir)
  ENDFOREACH(_file)
  
  INSTALL( CODE
      "EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_Fortran_MODULE_DIRECTORY}\" \"${CMAKE_INSTALL_PREFIX}/include/selalib\")")
   
   INSTALL(FILES ${SLL_MACROS}  DESTINATION include/selalib  )

ENDIF(BUILD_PACKAGE)
