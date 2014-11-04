##############################################################################
###      Sequential simulations
##############################################################################

ADD_LIBRARY( sll_simulations_sequential STATIC ../simulation_base_class.F90 )

TARGET_LINK_LIBRARIES( sll_simulations_sequential
                       sll_interpolators
                       sll_fields
                       sll_poisson_solvers
                       sll_reduction
                       sll_common_array_initializers
                     )

SET( INSTALL_DIR "${INSTALL_DIR}/simulations_sequential" )

ADD_SUBDIRECTORY( sim_VlasovPoisson1D_curvilinear_analytic )
ADD_SUBDIRECTORY( sim_VlasovPoisson1D_fourier_FEM          )

IF(HDF5_ENABLED)
  SET( SIM_LIBS sll_advection )
  ADD_SUBDIRECTORY( sim_2d_analytic_field_cartesian )
  ADD_SUBDIRECTORY( sim_2d_analytic_field_curvilinear )
  ADD_SUBDIRECTORY( sim_2d_guiding_center_cartesian )
  ADD_SUBDIRECTORY( sim_2d_guiding_center_polar )
ENDIF(HDF5_ENABLED)

