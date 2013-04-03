
SET(SLL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

SET(SLL_MACROS
	${CMAKE_CURRENT_SOURCE_DIR}/selalib.h
	${SLL_DIR}/memory/sll_memory.h
	${SLL_DIR}/precision/sll_working_precision.h
	${SLL_DIR}/assert/sll_assert.h
	${SLL_DIR}/constants/sll_constants.h
	${SLL_DIR}/utilities/sll_utilities.h
	${SLL_DIR}/maxwell_solvers/sll_maxwell_solvers.h
	${SLL_DIR}/poisson_solvers/sll_poisson_solvers.h
	)

SET(SLL_SOURCE_FILES
    ${SLL_DIR}/memory/sll_strings.F90
    ${SLL_DIR}/memory/sll_memory.F90
    ${SLL_DIR}/precision/sll_working_precision.F90
    ${SLL_DIR}/assert/sll_assert.F90
    ${SLL_DIR}/constants/sll_constants.F90
    ${SLL_DIR}/utilities/sll_utilities.F90
    ${SLL_DIR}/low_level_file_io/sll_ascii_io.F90
    ${SLL_DIR}/low_level_file_io/sll_binary_io.F90
    ${SLL_DIR}/low_level_file_io/sll_gnuplot.F90
    ${SLL_DIR}/low_level_file_io/sll_hdf5_io.F90
    ${SLL_DIR}/low_level_file_io/sll_xdmf.F90
    ${SLL_DIR}/low_level_file_io/sll_ascii_io.F90
    ${SLL_DIR}/low_level_file_io/sll_xml_io.F90
    ${SLL_DIR}/low_level_file_io/sll_io.F90
    ${SLL_DIR}/visu_pic/sll_visu_pic.F90
    ${SLL_DIR}/timer/sll_timer.F90
    ${SLL_DIR}/electric_field_accumulators/electric_field_2d.F90
    ${SLL_DIR}/splines/sll_splines.F90
    ${SLL_DIR}/splines/sll_cubic_splines.F90
    ${SLL_DIR}/splines/sll_quintic_splines.F90
    ${SLL_DIR}/splines/sll_odd_degree_splines.F90
    ${SLL_DIR}/splines/sll_arbitrary_degree_splines.F90
    ${SLL_DIR}/interpolators/sll_interpolators.F90
    ${SLL_DIR}/interpolators/sll_interpolators_1d_base.F90
    ${SLL_DIR}/interpolators/sll_odd_degree_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_quintic_spline_interpolator_1d.F90
    ${SLL_DIR}/interpolators/sll_cubic_spline_interpolator_2d.F90
    ${SLL_DIR}/interpolators/sll_interpolators_2d_base.F90
    ${SLL_DIR}/interpolators/sll_cubic_spline_interpolator_1d.F90
    ${SLL_DIR}/numerical_utilities/tridiagonal.F90
    ${SLL_DIR}/numerical_utilities/sll_toeplitz_penta_diagonal.F90
    ${SLL_DIR}/numerical_utilities/lagrange_interpolation.F90
    ${SLL_DIR}/maxwell_solvers/sll_maxwell_solvers.F90
    ${SLL_DIR}/maxwell_solvers/sll_maxwell_2d_fdtd.F90
    ${SLL_DIR}/maxwell_solvers/sll_maxwell_2d_pstd.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_solvers.F90
    ${SLL_DIR}/poisson_solvers/fft_module.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_1d_periodic.F90
    ${SLL_DIR}/poisson_solvers/sll_poisson_2d_periodic.F90)

ADD_LIBRARY(selalib STATIC ${SLL_SOURCE_FILES} selalib.F90)

SET(SLL_MODULES
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/selalib.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_memory.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_working_precision.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_assert.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_utilities.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_constants.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_maxwell_solvers.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_poisson_solvers.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_interpolators.mod
	 ${CMAKE_Fortran_MODULE_DIRECTORY}/sll_splines.mod
	)

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/selalib.h ${SLL_INCLUDE_FILES}
        DESTINATION include)
INSTALL(TARGETS selalib ARCHIVE DESTINATION lib)
INSTALL(FILES ${SLL_MACROS} DESTINATION include)
INSTALL(FILES ${SLL_MODULES} DESTINATION include)

IF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)

   ADD_LIBRARY(selalib-mpi STATIC ${SLL_SOURCE_FILES}
      ${SLL_DIR}/collective/sll_collective.F90
      ${SLL_DIR}/remap/sll_remap.F90
      ${SLL_DIR}/low_level_file_io/sll_hdf5_io_parallel.F90
      ${SLL_DIR}/low_level_file_io/sll_xdmf_parallel.F90 )

   INSTALL(TARGETS selalib-mpi ARCHIVE DESTINATION lib)
   FILE(COPY selalib-mpi.h DESTINATION ${CMAKE_BINARY_DIR}/include)
   INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/selalib-mpi.h   DESTINATION include)
   #INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/selalib-mpi.mod DESTINATION include)

ENDIF(HDF5_ENABLE_PARALLEL AND HDF5_PARALLEL_ENABLED)


## set the output directory for library and include
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${sll_BINARY_DIR}/lib)
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
#
## copy the main header files in the include directory 
#file(COPY
#       ${CMAKE_SOURCE_DIR}/selalib_level/selalib_low_level.h
#       ${CMAKE_SOURCE_DIR}/selalib_level/selalib_mid_level.h
#     DESTINATION
#       ${CMAKE_BINARY_DIR}/include/
#    )
#
#
##################################################################
## DEFINE LOW LEVEL STUFF
##################################################################
#set(LOW_LEVEL_INCLUDE_FILES
#	${CMAKE_SOURCE_DIR}/memory/sll_memory.h
#	${CMAKE_SOURCE_DIR}/precision/sll_working_precision.h
#	${CMAKE_SOURCE_DIR}/assert/sll_assert.h
#	${CMAKE_SOURCE_DIR}/utilities/sll_utilities.h
#	)
#
#set(LOW_LEVEL_SOURCE_FILES
#	${CMAKE_SOURCE_DIR}/memory/sll_strings.F90
#	${CMAKE_SOURCE_DIR}/memory/sll_memory.F90
#	${CMAKE_SOURCE_DIR}/precision/sll_working_precision.F90
#	${CMAKE_SOURCE_DIR}/assert/sll_assert.F90
#	${CMAKE_SOURCE_DIR}/constants/sll_constants.F90
#	${CMAKE_SOURCE_DIR}/utilities/utilities.F90
#	${CMAKE_SOURCE_DIR}/low_level_file_io/sll_ascii_io.F90
#	${CMAKE_SOURCE_DIR}/low_level_file_io/sll_binary_io.F90
#	${CMAKE_SOURCE_DIR}/low_level_file_io/sll_gnuplot.F90
#	${CMAKE_SOURCE_DIR}/low_level_file_io/sll_hdf5_io.F90
#	${CMAKE_SOURCE_DIR}/low_level_file_io/sll_io.F90
#	${CMAKE_SOURCE_DIR}/low_level_file_io/sll_xdmf.F90
#	${CMAKE_SOURCE_DIR}/low_level_file_io/sll_xml_io.F90
#	${CMAKE_SOURCE_DIR}/timer/sll_timer.F90
#	)
#
## copy all low level include files in the include directory
#file(COPY ${LOW_LEVEL_INCLUDE_FILES} DESTINATION ${CMAKE_BINARY_DIR}/include/)
#
#add_library(selalib STATIC ${LOW_LEVEL_SOURCE_FILES})
#
#set(HIGH_LEVEL_INCLUDE_FILES
#	${CMAKE_SOURCE_DIR}/remap/sll_remap.h
#	)
#
## Attention always use selalib fft, need to be fix!
## There are no poisson solver because the dependencies aren't obvious. 
#set(HIGH_LEVEL_SOURCE_FILES
#	${CMAKE_SOURCE_DIR}/numerical_utilities/tridiagonal.F90
#	${CMAKE_SOURCE_DIR}/non_linear_solvers/newton_raphson.F90
#	${CMAKE_SOURCE_DIR}/splines/util_constants.F90
#	${CMAKE_SOURCE_DIR}/splines/sll_splines.F90
#	${CMAKE_SOURCE_DIR}/splines/cubic_nonuniform_splines.F90
#	${CMAKE_SOURCE_DIR}/splines/sll_arbitrary_degree_splines.F90
#	${CMAKE_SOURCE_DIR}/integration/gauss_legendre.F90
#	${CMAKE_SOURCE_DIR}/fft/fft_utils.F90
#	${CMAKE_SOURCE_DIR}/fft/sll_fft_selalib_implementation.F90
#	${CMAKE_SOURCE_DIR}/collective/sll_collective.F90
#	${CMAKE_SOURCE_DIR}/remap/sll_remap.F90
#	${CMAKE_SOURCE_DIR}/WENO/WENO_interp.F90
#	${CMAKE_SOURCE_DIR}/WENO/WENO_recon.F90
#	)
#
## copy all law level include files in the include directory
#file(COPY ${HIGH_LEVEL_INCLUDE_FILES} DESTINATION ${CMAKE_BINARY_DIR}/include/)
#
####add_library(selalib_hl STATIC ${HIGH_LEVEL_SOURCE_FILES} ${LOW_LEVEL_SOURCE_FILES})
