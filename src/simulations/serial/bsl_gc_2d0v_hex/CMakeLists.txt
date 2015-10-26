IF (HDF5_ENABLED)
set(    SIM_NAME bsl_gc_2d0v_hex              )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} ${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} 
  sll_sim_serial
  sll_splines
  sll_meshes
  sll_hermite_interpolation
  hex_poisson
  sll_utilities)

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

# Ctest
# IF(BUILD_TESTING)
#   FILE( COPY ${INPUT_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#   ADD_TEST( NAME ${TEST_NAME1} COMMAND ${EXECUTABLE1} )
#   SET_TESTS_PROPERTIES(${TEST_NAME1} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
# ENDIF(BUILD_TESTING)

ENDIF (HDF5_ENABLED)





#set( SIM_NAME1   2d_guiding_center_hexagonal_splines )
#ADD_EXECUTABLE( ${EXECUTABLE1} simulation_${SIM_NAME1}.F90)
#TARGET_LINK_LIBRARIES( ${EXECUTABLE1} sll_sim_serial ${SIM_LIBS}
#                                                                sll_splines
#                                                                sll_meshes
#                                                                hex_poisson
#                                                                sll_utilities
#                                                                sll_hermite_interpolation)

#set( SIM_NAME2   2d_guiding_center_hexagonal_hermite )
#ADD_EXECUTABLE( ${EXECUTABLE2} simulation_${SIM_NAME2}.F90)
#TARGET_LINK_LIBRARIES( ${EXECUTABLE2} sll_sim_serial ${SIM_LIBS}
#                                                    sll_hermite_interpolation
#                                                    sll_meshes
#                                                    hex_poisson
#                                                    sll_utilities)

#set( SIM_NAME3   2d_rotation_hexagonal )
#ADD_EXECUTABLE( ${EXECUTABLE3} simulation_${SIM_NAME3}.F90)
#TARGET_LINK_LIBRARIES(
#  ${EXECUTABLE3}
#  sll_sim_serial
#  ${SIM_LIBS}
#  sll_hermite_interpolation
#  sll_splines
#  sll_meshes)

#set( SIM_NAME4   2d_rotation_polar )
#ADD_EXECUTABLE( ${EXECUTABLE4} simulation_${SIM_NAME4}.F90)
#TARGET_LINK_LIBRARIES(
#  ${EXECUTABLE4}
#  sll_sim_serial
#  ${SIM_LIBS}
#  sll_hermite_interpolation
#  sll_splines
#  sll_meshes)

#set( SIM_NAME5   2d_guiding_center_hexagonal )
#ADD_EXECUTABLE( ${EXECUTABLE5} simulation_${SIM_NAME5}.F90)
#TARGET_LINK_LIBRARIES( ${EXECUTABLE5} sll_sim_serial ${SIM_LIBS}
#                                                                sll_splines
#                                                                sll_meshes
#                                                                sll_hermite_interpolation
#                                                                hex_poisson
#                                                                sll_utilities)

