set(    SIM_NAME bsl_vp_1d1v_cart_serial      )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} ${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} 
  sll_interpolators
  sll_fft
  sll_distribution_function
  sll_file_io
  sll_utilities
  sll_poisson_solvers
  dfftpack )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )
