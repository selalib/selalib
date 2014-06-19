SET(SLL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

SET(SLL_MACROS
    ${CMAKE_CURRENT_SOURCE_DIR}/selalib.h
    ${CMAKE_CURRENT_SOURCE_DIR}/selalib-mpi.h
    ${SLL_DIR}/memory/sll_memory.h
    ${SLL_DIR}/precision/sll_working_precision.h
    ${SLL_DIR}/assert/sll_assert.h
    ${SLL_DIR}/constants/sll_constants.h
    ${SLL_DIR}/utilities/sll_utilities.h
    ${SLL_DIR}/maxwell_solvers/sll_maxwell_solvers_macros.h
    ${SLL_DIR}/maxwell_solvers/sll_maxwell_solvers.h
    ${SLL_DIR}/poisson_solvers/sll_poisson_solvers.h
    ${SLL_DIR}/interpolators/sll_interpolators.h
    ${SLL_DIR}/splines/sll_splines.h
    ${SLL_DIR}/file_io/sll_file_io.h
    ${SLL_DIR}/fft/sll_fft.h
    ${SLL_DIR}/fft/sll_fftw.h
    ${SLL_DIR}/integration/sll_integration.h
    ${SLL_DIR}/fields/sll_field_1d.h
    ${SLL_DIR}/fields/sll_field_2d.h
    ${SLL_DIR}/boundary_condition_descriptors/sll_boundary_condition_descriptors.h
    ${SLL_DIR}/logical_meshes/sll_logical_meshes.h
    ${SLL_DIR}/fields/sll_fields.h
    ${SLL_DIR}/mapped_meshes/sll_mapped_meshes.h
    ${SLL_DIR}/coordinate_transformations/sll_coordinate_transformations.h
    ${SLL_DIR}/pastix/sll_pastix.h
    ${SLL_DIR}/linear_solvers/sll_linear_solvers.h
    ${SLL_DIR}/pic_accumulators/sll_accumulators.h
    ${SLL_DIR}/pic_particle_types/particle_representation.h
)

FILE(COPY ${SLL_MACROS} DESTINATION ${CMAKE_BINARY_DIR}/include)

SET(SLL_SOURCE_FILES
    ${SLL_DIR}/memory/sll_strings.F90
    ${SLL_DIR}/memory/sll_memory.F90
    ${SLL_DIR}/precision/sll_working_precision.F90
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
    ${SLL_DIR}/constants/sll_constants.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_explicit_euler_conservative.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_trapezoid.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_trapezoid_conservative.F90
    ${SLL_DIR}/characteristics/sll_characteristics_1d_base.F90
    ${SLL_DIR}/characteristics/sll_characteristics_2d_verlet.F90
    ${SLL_DIR}/characteristics/sll_characteristics_2d_explicit_euler.F90
    ${SLL_DIR}/characteristics/sll_characteristics_2d_explicit_euler_conservative.F90
    ${SLL_DIR}/characteristics/sll_characteristics_2d_base.F90
    ${SLL_DIR}/utilities/sll_utilities.F90
    ${SLL_DIR}/utilities/tridiagonal.F90
    ${SLL_DIR}/utilities/sll_toeplitz_penta_diagonal.F90
    ${SLL_DIR}/lagrange_interpolation/sll_lagrange_interpolation.F90
    ${SLL_DIR}/file_io/sll_file_io.F90
    ${SLL_DIR}/file_io/sll_binary_io.F90
    ${SLL_DIR}/file_io/sll_gnuplot.F90
    ${SLL_DIR}/file_io/sll_hdf5_io.F90
    ${SLL_DIR}/file_io/sll_xdmf.F90
    ${SLL_DIR}/file_io/sll_ascii_io.F90
    ${SLL_DIR}/file_io/sll_xml_io.F90
    ${SLL_DIR}/file_io/sll_plotmtv.F90
    ${SLL_DIR}/visu_pic/sll_visu_pic.F90
    ${SLL_DIR}/timer/sll_timer.F90
    ${SLL_DIR}/field_accumulators/sll_field_accumulators_2d.F90
    ${SLL_DIR}/fields/sll_scalar_field_2d_base.F90
    ${SLL_DIR}/fields/sll_scalar_field_2d_alternative.F90
    ${SLL_DIR}/general_coordinate_elliptic_solvers/sll_general_coordinate_elliptic_solver.F90
    ${SLL_DIR}/general_coordinate_elliptic_solvers/sll_knots.F90
    ${SLL_DIR}/general_coordinate_elliptic_solvers/connectivity_qn.f90
    ${SLL_DIR}/splines/sll_cubic_splines.F90
    ${SLL_DIR}/splines/sll_quintic_splines.F90
    ${SLL_DIR}/splines/sll_odd_degree_splines.F90
    ${SLL_DIR}/splines/sll_arbitrary_degree_splines.F90
    ${SLL_DIR}/splines/cubic_nonuniform_splines.F90
    ${SLL_DIR}/deboor_splines/dvalue1d.f90
    ${SLL_DIR}/deboor_splines/bvalue2d.f90
    ${SLL_DIR}/deboor_splines/dvalue2d.f90
    ${SLL_DIR}/deboor_splines/bsplvb.f90
    ${SLL_DIR}/deboor_splines/bvalue.f90
    ${SLL_DIR}/deboor_splines/spli1d_per.f90
    ${SLL_DIR}/deboor_splines/spli1d_dir.f90
    ${SLL_DIR}/deboor_splines/spli2d_per.f90
    ${SLL_DIR}/deboor_splines/spli2d.f90
    ${SLL_DIR}/deboor_splines/spli2d_custom.f90
    ${SLL_DIR}/deboor_splines/spli2d_perper.f90
    ${SLL_DIR}/deboor_splines/banfac.f90  
    ${SLL_DIR}/deboor_splines/interv.f90  
    ${SLL_DIR}/deboor_splines/splint.f90
    ${SLL_DIR}/deboor_splines/banslv.f90  
    ${SLL_DIR}/deboor_splines/bsplvd.f90
    ${SLL_DIR}/deboor_splines/knots.f90
    ${SLL_DIR}/deboor_splines/splopt.f90
    ${SLL_DIR}/deboor_splines/sll_module_deboor_splines_1d.F90
    ${SLL_DIR}/deboor_splines/sll_module_deboor_splines_2d.F90
    ${SLL_DIR}/interpolators/sll_interpolators_1d_base.F90
    ${SLL_DIR}/interpolators/sll_odd_degree_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_quintic_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_cubic_spline_interpolator_2d.F90
    ${SLL_DIR}/interpolators/sll_interpolators_2d_base.F90
    ${SLL_DIR}/interpolators/sll_cubic_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_arbitrary_degree_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_arbitrary_degree_spline_interpolator_2d.F90
    ${SLL_DIR}/maxwell_solvers/sll_maxwell_solvers_base.F90
    ${SLL_DIR}/maxwell_solvers/sll_maxwell_2d_fdtd.F90
    ${SLL_DIR}/periodic_interpolation/periodic_interp.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_solvers.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_1d_periodic.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_2d_polar.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_2d_periodic_fftpack.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_1d_base.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_1d_periodic_solver.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_1d_polar_solver.F90
    ${SLL_DIR}/fft/sll_fft_utils.F90
    ${SLL_DIR}/integration/gauss_legendre.F90
    ${SLL_DIR}/integration/gauss_lobatto.F90
    ${SLL_DIR}/simulation/sll_polar_operators.F90
    ${SLL_DIR}/simulation/sll_advection_polar.F90
#    ${SLL_DIR}/sparse_matrix_interface/basic_cg.F90
    ${SLL_DIR}/sparse_matrix_interface/gauss_seidel.F90
    ${SLL_DIR}/sparse_matrix_interface/getdia.F90
    ${SLL_DIR}/sparse_matrix_interface/getl.F90
    ${SLL_DIR}/sparse_matrix_interface/getu.F90
    ${SLL_DIR}/sparse_matrix_interface/lsol.F90
    ${SLL_DIR}/sparse_matrix_interface/remove_zeros.F90
    ${SLL_DIR}/sparse_matrix_interface/sll_sparse_matrix_module.F90
#    ${SLL_DIR}/sparse_matrix_interface/sll_sparse_matrix_module_basic.F90
    ${SLL_DIR}/sparse_matrix_interface/spm_basic.F90
    ${SLL_DIR}/sparse_matrix_interface/spm_def.F90
#    ${SLL_DIR}/sparse_matrix_interface/spm_mkl.F90
    ${SLL_DIR}/sparse_matrix_interface/spm_murge.F90
#    ${SLL_DIR}/sparse_matrix_interface/spm_template.F90
    ${SLL_DIR}/sparse_matrix_interface/spm_tools.F90
    ${SLL_DIR}/time_solvers/sll_time_splitting_coeff.F90
    ${SLL_DIR}/boundary_condition_descriptors/sll_boundary_condition_descriptors.F90
    ${SLL_DIR}/logical_meshes/sll_logical_meshes.F90
    ${SLL_DIR}/fields/scalar_field_initializer_base.F90
#    ${SLL_DIR}/fields/sll_scalar_field_2d_base.F90
    ${SLL_DIR}/fields/sll_scalar_field_2d.F90
    ${SLL_DIR}/mapped_meshes/sll_mapped_meshes_2d_base.F90
    ${SLL_DIR}/mapped_meshes/sll_mapped_meshes_2d.F90
    ${SLL_DIR}/coordinate_transformations/sll_common_coordinate_transformations.F90	
    ${SLL_DIR}/coordinate_transformations/sll_coordinate_transformation_2d_base.F90	
    ${SLL_DIR}/coordinate_transformations/sll_coordinate_transformations_2d.F90
    ${SLL_DIR}/coordinate_transformations/sll_coordinate_transformations_2d_nurbs.F90
    ${SLL_DIR}/mesh_calculus/sll_mesh_calculus_module.F90
 )

IF(Fortran_COMPILER MATCHES "IBM")
ELSE()
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/utilities/sol.f
    ${SLL_DIR}/utilities/muls.f)
ENDIF()

IF(${FFT_LIB} MATCHES "FFTPACK")

  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
      ${SLL_DIR}/fft/sll_fft_mod_through_fftpack.F90

      )
ELSEIF(${FFT_LIB} MATCHES "FFTW")

  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
      ${SLL_DIR}/fft/sll_fft_mod_through_fftw.F90
      ${SLL_DIR}/fft/fftw3.F90 )

ELSE()

  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
      ${SLL_DIR}/fft/sll_fft_selalib_implementation.F90)

ENDIF()

# Got the impression that this maxwell solver is made to work with FFTW only,
# thus I'm including the following. Intent needs to be made more explicit if 
# this is not the case. ECG
IF((${FFT_LIB} MATCHES "FFTW") AND FFTW_ENABLED)
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/maxwell_solvers/sll_maxwell_2d_pstd.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_2d_periodic_fftw.F90
  )
ENDIF()

IF(MUDPACK_ENABLED)

  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES} 
    ${SLL_DIR}/poisson_solvers/sll_mudpack_base.F90
    ${SLL_DIR}/poisson_solvers/sll_mudpack_cartesian.F90
  )

ENDIF(MUDPACK_ENABLED)


IF(SUITESPARSE_ENABLED)
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
        ${SLL_DIR}/sparse_matrix_interface/spm_umfpack.F90
	${SLL_DIR}/sparse_matrix_interface/sll_sparse_matrix_umfpack_module.F90
  )
ELSE()
  SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
      ${SLL_DIR}/sparse_matrix_interface/basic_cg.F90)

ENDIF(SUITESPARSE_ENABLED)


IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

   SET(SLL_SOURCE_FILES ${SLL_SOURCE_FILES}
       ${SLL_DIR}/collective/sll_collective.F90
       ${SLL_DIR}/remap/sll_remap.F90
       ${SLL_DIR}/point_to_point_communications/sll_point_to_point_comms.F90
       ${SLL_DIR}/file_io/sll_hdf5_io_parallel.F90
       ${SLL_DIR}/file_io/sll_xdmf_parallel.F90
       ${SLL_DIR}/file_io/sll_gnuplot_parallel.F90
       ${SLL_DIR}/poisson_solvers/sll_poisson_2d_periodic_cartesian_par.F90
       ${SLL_DIR}/poisson_solvers/sll_poisson_polar_parallel.F90
       ${SLL_DIR}/simulation/simulation_base_class.F90
       ${SLL_DIR}/parallel_array_initializers/sll_common_array_initializers_module.F90
       ${SLL_DIR}/parallel_array_initializers/sll_parallel_array_initializer_module.F90
       ${SLL_DIR}/parallel_utilities/buffer_loader_utilities.F90
       ${SLL_DIR}/simulation/simulation_4d_qns_general.F90
       ${SLL_DIR}/simulation/simulation_2d_vlasov_poisson_cartesian.F90
 )

ENDIF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

IF(SLL_BUILD_PACKAGE)

SET(BUILD_SHARED_LIBRARY  OFF CACHE BOOL "build shared library")
#SET(CMAKE_SHARED_LINKER_FLAGS "-fpic")
IF(BUILD_SHARED_LIBRARY)
   ADD_LIBRARY(selalib SHARED ${SLL_SOURCE_FILES})
   TARGET_LINK_LIBRARIES(selalib dfftpack mudpack ${HDF5_LIBRARIES} ${LAPACK_LIBRARIES})
ELSE()
   ADD_LIBRARY(selalib STATIC ${SLL_SOURCE_FILES})
   INSTALL(TARGETS selalib ARCHIVE DESTINATION lib COMPONENT libraries)
ENDIF()

SET(DIST_NAME "selalib-0.5.0")

ADD_CUSTOM_TARGET(release 
COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/${DIST_NAME}/include
COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/macports/README ${CMAKE_BINARY_DIR}/${DIST_NAME}/
COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/macports/CMakeLists.txt ${CMAKE_BINARY_DIR}/${DIST_NAME}/
COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/macports/Doxyfile.in ${CMAKE_BINARY_DIR}/${DIST_NAME}/Doxyfile
COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/../documentation.f90 ${CMAKE_BINARY_DIR}/${DIST_NAME}/
COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/../../documentation/sphinx/_static/logo2.png ${CMAKE_BINARY_DIR}/${DIST_NAME}/
COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../fftpack ${CMAKE_BINARY_DIR}/${DIST_NAME}/fftpack
COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/${DIST_NAME}/examples
COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/${DIST_NAME}/CMakeFiles
COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/${DIST_NAME}/CMakeCache.txt
COMMAND ${CMAKE_COMMAND} -E tar zcf ${DIST_NAME}.tar.gz ${CMAKE_BINARY_DIR}/${DIST_NAME}
COMMAND ${CMAKE_COMMAND} -E md5sum ${DIST_NAME}.tar.gz
COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/${DIST_NAME}
WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
COMMENT "Generating Release for macports" VERBATIM)

ADD_CUSTOM_TARGET(make_source_dir COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/${DIST_NAME}/src)
ADD_DEPENDENCIES(release make_source_dir)

FOREACH(_file ${SLL_SOURCE_FILES})
   get_filename_component(_target ${_file} NAME_WE)
   ADD_CUSTOM_TARGET(${_target}_f90 ${CMAKE_COMMAND} -E copy ${_file} ${CMAKE_BINARY_DIR}/${DIST_NAME}/src)
   ADD_DEPENDENCIES(release ${_target}_f90)
   ADD_DEPENDENCIES(${_target}_f90 make_source_dir)
ENDFOREACH(_file)

INSTALL( CODE
   "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_Fortran_MODULE_DIRECTORY}\" \"${CMAKE_INSTALL_PREFIX}/include/selalib\")")

INSTALL(FILES ${SLL_MACROS}  DESTINATION include/selalib  )

ENDIF(SLL_BUILD_PACKAGE)
