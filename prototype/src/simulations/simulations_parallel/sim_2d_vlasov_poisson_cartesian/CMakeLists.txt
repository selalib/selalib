set(   TEST_NAME 2d_vp_cartesian               )
set(    SIM_NAME 2d_vlasov_poisson_cartesian   )
set(  EXECUTABLE "test_${TEST_NAME}"           )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} simulation_${SIM_NAME}.F90
                              unit_test_${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} sll_simulations_parallel )

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Quick fix for ctest
FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_BINARY_DIR} )

#==============================================================================
# Postprocessing code for KEEN
#==============================================================================
IF(Fortran_COMPILER MATCHES IBM)
   SET_SOURCE_FILES_PROPERTIES( VPpostprocessing_drive_KEEN.f 
                                PROPERTIES COMPILE_FLAGS -qfixed )
ENDIF(Fortran_COMPILER MATCHES IBM)

ADD_EXECUTABLE( VPpostprocessing_drive_KEEN VPpostprocessing_drive_KEEN.f )
TARGET_LINK_LIBRARIES( VPpostprocessing_drive_KEEN dfftpack )

INSTALL( TARGETS VPpostprocessing_drive_KEEN DESTINATION ${INSTALL_DIR} )
INSTALL( FILES vp_plots.py DESTINATION ${INSTALL_DIR} )

