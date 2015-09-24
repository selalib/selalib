##############################################################################
###      Sequential simulations
##############################################################################

ADD_LIBRARY( sll_simulations_sequential STATIC ../simulation_base_class.F90 )

TARGET_LINK_LIBRARIES( sll_simulations_sequential
                       sll_interpolators
                       sll_fields
                       sll_poisson_solvers
                       sll_reduction
                       sll_operator_splitting
                       sll_common_array_initializers
                       sll_gyroaverage
                     )

SET( INSTALL_DIR "${INSTALL_DIR}/simulations_sequential" )

ADD_SUBDIRECTORY( sim_VlasovPoisson1D_curvilinear_analytic  )
ADD_SUBDIRECTORY( sim_VlasovPoisson1D_fourier_FEM           )
ADD_SUBDIRECTORY( sim_VlasovPoisson1D_cartesian_non_uniform )
ADD_SUBDIRECTORY( sim_VlasovPoisson1D_various_programs      )
ADD_SUBDIRECTORY( sim_guiding_center_2D_generalized_coords  )
ADD_SUBDIRECTORY( sim_polar_SL                              )
ADD_SUBDIRECTORY( sim_curvilinear_SL                        )
ADD_SUBDIRECTORY( sim_CG_polar                              )
ADD_SUBDIRECTORY( sim_CG_curvilinear_2D                     )
ADD_SUBDIRECTORY( sim_poisson_polar_test                    )
ADD_SUBDIRECTORY( sim_2d_guiding_center_hexagonal           )

IF(HDF5_ENABLED)
  SET( SIM_LIBS sll_advection sll_parallel_array_utilities )
  ADD_SUBDIRECTORY( sim_2d_analytic_field_cartesian )
  ADD_SUBDIRECTORY( sim_2d_analytic_field_curvilinear )
  ADD_SUBDIRECTORY( sim_2d_guiding_center_cartesian )
  ADD_SUBDIRECTORY( sim_2d_guiding_center_polar )
  ADD_SUBDIRECTORY( sim_2d_guiding_center_polar_one_mu )
  ADD_SUBDIRECTORY( sim_guiding_center_2D_generalized_coords_mudpack  )
ENDIF(HDF5_ENABLED)

ADD_SUBDIRECTORY(bsl_examples)
ADD_SUBDIRECTORY(csl_examples)
ADD_SUBDIRECTORY(fsl_examples)
