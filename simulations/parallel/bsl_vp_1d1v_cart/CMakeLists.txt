set(    SIM_NAME bsl_vp_1d1v_cart             )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

add_executable( ${EXECUTABLE} sll_m_sim_${SIM_NAME}.F90 ${TEST_NAME}.F90 )
target_include_directories(${EXECUTABLE} BEFORE PRIVATE ${HDF5_INCLUDE_DIR_FORTRAN})
target_link_libraries( ${EXECUTABLE} ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS} )

install( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
install( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

if(BUILD_TESTING)
  file( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  set(  ARGS ${CMAKE_CURRENT_BINARY_DIR}/vpsim2d_cartesian_input )
  add_mpi_test( ${TEST_NAME} ${EXECUTABLE} ${PROCS} ${ARGS} )
  set_tests_properties( ${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
endif(BUILD_TESTING)

#==============================================================================
# Postprocessing code for KEEN
#==============================================================================
if(Fortran_COMPILER MATCHES IBM)
   set_source_files_properties( VPpostprocessing_drive_KEEN.f 
                                PROPERTIES COMPILE_FLAGS -qfixed )
endif(Fortran_COMPILER MATCHES IBM)

add_executable( VPpostprocessing_drive_KEEN VPpostprocessing_drive_KEEN.f )
target_link_libraries( VPpostprocessing_drive_KEEN dfftpack )

install( TARGETS VPpostprocessing_drive_KEEN DESTINATION ${INSTALL_DIR} )
install( FILES vp_plots.py DESTINATION ${INSTALL_DIR} )
