SET(SLL_DIR ${CMAKE_SOURCE_DIR}/src)

FILE(GLOB_RECURSE SLL_MACROS ${PROJECT_SOURCE_DIR}/*.h )

SET(SLL_SOURCE_FILES
    ${SLL_DIR}/low_level_utilities/memory/sll_m_memory.F90
    ${SLL_DIR}/low_level_utilities/working_precision/sll_m_working_precision.F90
    ${SLL_DIR}/low_level_utilities/assert/sll_m_assert.F90
    ${SLL_DIR}/low_level_utilities/timer/sll_m_timer.F90
    ${SLL_DIR}/low_level_utilities/errors/sll_m_errors.F90
    ${SLL_DIR}/low_level_utilities/utilities/sll_m_utilities.F90
    ${SLL_DIR}/quadrature/sll_m_primitives.F90
    ${SLL_DIR}/low_level_utilities/utilities/sll_m_fornberg.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_1d_periodic.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_1d_non_uniform_cubic_splines.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_1d_BSL.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_1d_CSL.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_1d_PSM.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_1d_base.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_2d_bsl.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_2d_tensor_product.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_2d_CSL.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_2d_base.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_1d_ampere.F90
    ${SLL_DIR}/semi_lagrangian/advection/sll_m_advection_1d_spectral.F90
    ${SLL_DIR}/low_level_utilities/constants/sll_m_constants.F90
    ${SLL_DIR}/time_integration/splitting_methods/sll_m_time_splitting.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_1d_explicit_euler.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_1d_explicit_euler_conservative.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_1d_trapezoid.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_1d_trapezoid_conservative.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_1d_base.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_2d_verlet.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_2d_explicit_euler.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_2d_explicit_euler_conservative.F90
    ${SLL_DIR}/time_integration/characteristics/sll_m_characteristics_2d_base.F90
    ${SLL_DIR}/linear_solvers/sll_m_tridiagonal.F90
    ${SLL_DIR}/linear_solvers/sll_m_penta_diagonal.F90
    ${SLL_DIR}/../external/pppack/banfac.f90
    ${SLL_DIR}/../external/pppack/banslv.f90
    ${SLL_DIR}/interpolation/lagrange_interpolation/sll_m_lagrange_interpolation_1d.F90
    ${SLL_DIR}/io/file_io/sll_m_binary_io.F90
    ${SLL_DIR}/io/file_io/sll_m_gnuplot.F90
    ${SLL_DIR}/io/file_io/sll_m_hdf5_io_serial.F90
    ${SLL_DIR}/io/file_io/sll_m_xdmf.F90
    ${SLL_DIR}/io/file_io/sll_m_ascii_io.F90
    ${SLL_DIR}/io/file_io/sll_m_xml_io.F90
    ${SLL_DIR}/io/file_io/sll_m_plotmtv.F90
    ${SLL_DIR}/io/file_io/sll_m_tri_mesh_xmf.F90
    ${SLL_DIR}/particle_methods/pic_visu/sll_m_pic_visu.F90
    ${SLL_DIR}/data_structures/fields/sll_m_scalar_field_1d_base.F90
    ${SLL_DIR}/data_structures/fields/sll_m_scalar_field_1d.F90
    ${SLL_DIR}/data_structures/fields/sll_m_scalar_field_2d_base.F90
    ${SLL_DIR}/data_structures/fields/sll_m_scalar_field_2d.F90
    ${SLL_DIR}/field_solvers/general_coordinate_elliptic_solvers/sll_m_general_coordinate_elliptic_solver.F90
    ${SLL_DIR}/splines/sll_m_bsplines.F90
    ${SLL_DIR}/splines/sll_m_cubic_splines.F90
    ${SLL_DIR}/splines/sll_m_arbitrary_degree_splines.F90
    ${SLL_DIR}/splines/sll_m_cubic_non_uniform_splines.F90
    ${SLL_DIR}/splines/sll_m_deboor_splines_1d.F90
    ${SLL_DIR}/splines/sll_m_deboor_splines_2d.F90
    ${SLL_DIR}/interpolation/hermite_interpolation/sll_m_hermite_interpolation_2d.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_interpolators_1d_base.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_interpolators_2d_base.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_cubic_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_cubic_spline_interpolator_1d_nonuniform.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_cubic_spline_interpolator_2d.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_arbitrary_degree_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_arbitrary_degree_spline_interpolator_2d.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_periodic_interpolator_1d.F90
    ${SLL_DIR}/interpolation/interpolators/sll_m_lagrange_interpolator_1d.F90
    ${SLL_DIR}/field_solvers/maxwell_solvers/sll_m_maxwell_solvers_base.F90
    ${SLL_DIR}/field_solvers/maxwell_solvers/sll_m_maxwell_2d_fdtd.F90
    ${SLL_DIR}/interpolation/periodic_interpolation/sll_m_periodic_interp.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_1d_periodic.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_2d_polar.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_2d_periodic_fftpack.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_1d_base.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_1d_periodic_solver.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_1d_polar_solver.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_2d_polar_wrapper.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_2d_base.F90
    ${SLL_DIR}/interfaces/fft/sll_m_fft_utils.F90
    ${SLL_DIR}/interfaces/mpi/sll_mpi.F90
    ${SLL_DIR}/quadrature/sll_m_rectangle_integration.F90
    ${SLL_DIR}/quadrature/sll_m_trapz_integration.F90
    ${SLL_DIR}/quadrature/sll_m_gauss_legendre_integration.F90
    ${SLL_DIR}/quadrature/sll_m_gauss_lobatto_integration.F90
    #${SLL_DIR}/simulations/serial/sim_CG_polar/sll_polar_operators.F90
    #${SLL_DIR}/simulations/serial/sim_CG_polar/sll_advection_polar.F90
    ${SLL_DIR}/interfaces/sparse_matrix_manager/sll_m_qsort_partition.F90
    ${SLL_DIR}/interfaces/sparse_matrix_manager/sll_m_sparse_matrix.F90
    ${SLL_DIR}/time_integration/splitting_methods/sll_m_time_splitting_coeff.F90
    ${SLL_DIR}/data_structures/boundary_condition_descriptors/sll_m_boundary_condition_descriptors.F90
    ${SLL_DIR}/mesh/meshes/sll_m_meshes_base.F90
    ${SLL_DIR}/mesh/meshes/sll_m_cartesian_meshes.F90
    ${SLL_DIR}/mesh/meshes/sll_m_triangular_meshes.F90
    ${SLL_DIR}/mesh/meshes/sll_m_hexagonal_meshes.F90
    ${SLL_DIR}/data_structures/fields/sll_m_parallel_field_initializers_base.F90
    ${SLL_DIR}/data_structures/distribution_function/sll_m_scalar_field_initializers_base.F90
    ${SLL_DIR}/data_structures/distribution_function/sll_m_scalar_field_2d_old.F90
    ${SLL_DIR}/mesh/coordinate_transformations/sll_m_common_coordinate_transformations.F90
    ${SLL_DIR}/mesh/coordinate_transformations/sll_m_coordinate_transformations_2d_nurbs.F90
    ${SLL_DIR}/mesh/coordinate_transformations/sll_m_coordinate_transformation_2d_base.F90	
    ${SLL_DIR}/mesh/coordinate_transformations/sll_m_coordinate_transformations_2d.F90
    ${SLL_DIR}/add_ons/multipatch/sll_m_cartesian_meshes_multipatch.F90
    ${SLL_DIR}/add_ons/multipatch/sll_m_coordinate_transformation_multipatch.F90
    ${SLL_DIR}/add_ons/multipatch/sll_m_scalar_field_2d_multipatch.F90
    ${SLL_DIR}/add_ons/multipatch/sll_m_general_coordinate_elliptic_solver_multipatch.F90
    ${SLL_DIR}/interfaces/sparse_matrix_manager/sll_m_sparse_matrix_mp.F90
    ${SLL_DIR}/mesh/mesh_calculus/sll_m_mesh_calculus_2d.F90
    ${SLL_DIR}/data_structures/distribution_function/sll_m_distribution_function_initializer_4d.F90
    ${SLL_DIR}/data_structures/distribution_function/sll_m_fdistribu4d_dk.F90
    ${SLL_DIR}/parallelization/reduction/sll_m_reduction.F90
    ${SLL_DIR}/particle_methods/random_deviate_generators/sll_m_gaussian.F90
    ${SLL_DIR}/particle_methods/random_deviate_generators/sll_m_hammersley.F90
    ${SLL_DIR}/particle_methods/random_deviate_generators/moment_matching/sll_m_moment_matching.F90
    ${SLL_DIR}/../external/burkardt/sobol/sll_m_sobol.F90
    ${SLL_DIR}/../external/burkardt/sobol/sobol.f90
    ${SLL_DIR}/../external/burkardt/prob/sll_m_prob.F90
    ${SLL_DIR}/../external/burkardt/prob/prob.f90
)

IF(F2008)
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/data_structures/vector_space/sll_m_vector_space_base.F90
    ${SLL_DIR}/data_structures/vector_space/sll_m_vector_space_real_arrays.F90
    ${SLL_DIR}/time_integration/ode_integrators/sll_m_ode_integrator_base.F90
    ${SLL_DIR}/time_integration/ode_integrators/sll_m_rk_explicit.F90)
ENDIF(F2008)

IF(Fortran_COMPILER MATCHES "IBM")
ELSE()
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/linear_solvers/sol.f
    ${SLL_DIR}/linear_solvers/muls.f)
ENDIF()

IF(${FFT_LIB} MATCHES "FFTPACK")
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} ${SLL_DIR}/interfaces/fft/sll_m_fft_fftpack.F90)
ELSEIF(${FFT_LIB} MATCHES "FFTW")
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} ${SLL_DIR}/interfaces/fft/sll_m_fft_fftw3.F90) 
ELSE()
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} ${SLL_DIR}/interfaces/fft/sll_m_fft_sllfft.F90)
ENDIF()

IF(FFTW_ENABLED)
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/interfaces/fft/sll_m_fftw3.F90
    ${SLL_DIR}/field_solvers/maxwell_solvers/sll_m_maxwell_2d_pstd.F90
    ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_poisson_2d_periodic_fftw.F90
  )
ENDIF()

SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
  ${SLL_DIR}/field_solvers/poisson_solvers/sll_m_mudpack.F90
)

IF(SUITESPARSE_ENABLED)
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
    ${SLL_DIR}/interfaces/sparse_matrix_manager/spm_umfpack.F90
    ${SLL_DIR}/interfaces/sparse_matrix_manager/sll_sparse_matrix_umfpack_module.F90
  )
ENDIF(SUITESPARSE_ENABLED)

IF(MPI_ENABLED)
 SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
  ${SLL_DIR}/parallelization/collective/sll_m_collective.F90
  ${SLL_DIR}/parallelization/remap/sll_m_remapper.F90
  ${SLL_DIR}/parallelization/point_to_point_communications/sll_m_point_to_point_comms.F90
  ${SLL_DIR}/field_solvers/poisson_solvers_parallel/sll_m_poisson_2d_periodic_cartesian_par.F90
  ${SLL_DIR}/field_solvers/poisson_solvers_parallel/sll_m_poisson_polar_parallel.F90
  ${SLL_DIR}/simulations/simulation_base_class.F90
  ${SLL_DIR}/initialization_tools/function_initialization/sll_m_common_array_initializers.F90
  ${SLL_DIR}/parallelization/parallel_array_utilities/sll_m_parallel_array_initializer.F90
  ${SLL_DIR}/parallelization/parallel_utilities/sll_m_buffer_loader_utilities.F90
  ${SLL_DIR}/io/file_io_parallel/sll_m_gnuplot_parallel.F90)
ENDIF(MPI_ENABLED)

IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

 SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
  ${SLL_DIR}/io/file_io_parallel/sll_m_hdf5_io_parallel.F90
  ${SLL_DIR}/io/file_io_parallel/sll_m_xdmf_parallel.F90
  ${SLL_DIR}/add_ons/multipatch/sll_m_distribution_function_4d_multipatch.F90
  ${SLL_DIR}/simulations/parallel/bsl_vp_1d1v_cart/sll_m_sim_bsl_vp_1d1v_cart.F90
  ${SLL_DIR}/simulations/parallel/bsl_va_1d1v_cart_spectral/sll_m_sim_bsl_va_1d1v_cart_spectral.F90
  ${SLL_DIR}/simulations/parallel/bsl_vp_2d2v_cart_multipatch/sll_m_sim_bsl_vp_2d2v_cart_multipatch.F90
  ${SLL_DIR}/simulations/parallel/bsl_dk_3d1v_curv/sll_m_sim_bsl_dk_3d1v_curv.F90
  ${SLL_DIR}/simulations/parallel/bsl_dk_3d1v_polar/sll_m_sim_bsl_dk_3d1v_polar.F90)

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
          ${CMAKE_SOURCE_DIR}/src/documentation.f90 
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
          ${CMAKE_SOURCE_DIR}/external/fftpack 
          ${CMAKE_BINARY_DIR}/${DIST_NAME}/fftpack
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
          ${CMAKE_SOURCE_DIR}/external/mudpack 
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
