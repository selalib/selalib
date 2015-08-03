IF (HDF5_ENABLED)
# Whole simulation :
set( TEST_NAME 2d_gc_hex )
set( INSTALL_DIR "${INSTALL_DIR}/${TEST_NAME}" )
set( INPUT_FILE sim2d_gc_hex_input.nml )

# First test :
set( TEST_NAME1  2d_gc_hex_splines )
set( SIM_NAME1   2d_guiding_center_hexagonal_splines )
set( EXECUTABLE1 "test_${TEST_NAME1}" )

# Second test :
set( TEST_NAME2  2d_gc_hex_hermite )
set( SIM_NAME2   2d_guiding_center_hexagonal_hermite )
set( EXECUTABLE2 "test_${TEST_NAME2}" )

# Third test :
set( TEST_NAME3  2d_rot_hex_splines )
set( SIM_NAME3   2d_rotation_hexagonal_splines )
set( EXECUTABLE3 "test_${TEST_NAME3}" )

# Fourth test :
set( TEST_NAME4  2d_rot_hex )
set( SIM_NAME4   2d_rotation_hexagonal )
set( EXECUTABLE4 "test_${TEST_NAME4}" )

# Fifth test :
set( TEST_NAME5  2d_rot_polar )
set( SIM_NAME5   2d_rotation_polar )
set( EXECUTABLE5 "test_${TEST_NAME5}" )



# First test :
ADD_EXECUTABLE( ${EXECUTABLE1} simulation_${SIM_NAME1}.F90)
TARGET_LINK_LIBRARIES( ${EXECUTABLE1} sll_simulations_sequential ${SIM_LIBS}
                                                                sll_box_splines
                                                                sll_hex_meshes
                                                                hex_poisson
                                                                pivotbande
                                                                euler_2d_hex)

# Second test :
ADD_EXECUTABLE( ${EXECUTABLE2} simulation_${SIM_NAME2}.F90)
TARGET_LINK_LIBRARIES( ${EXECUTABLE2} sll_simulations_sequential ${SIM_LIBS}
                                                                hex_hermite
                                                                sll_hex_meshes
                                                                hex_poisson
                                                                pivotbande
                                                                euler_2d_hex)

# Third test :
ADD_EXECUTABLE( ${EXECUTABLE3} simulation_${SIM_NAME3}.F90)
TARGET_LINK_LIBRARIES(
  ${EXECUTABLE3} 
  sll_simulations_sequential
  ${SIM_LIBS}
  sll_box_splines
  sll_hex_meshes)

# Fourth test :
ADD_EXECUTABLE( ${EXECUTABLE4} simulation_${SIM_NAME4}.F90)
TARGET_LINK_LIBRARIES(
  ${EXECUTABLE4} 
  sll_simulations_sequential
  ${SIM_LIBS}
  hex_hermite
  sll_box_splines
  sll_hex_meshes)

# Fifth test :
ADD_EXECUTABLE( ${EXECUTABLE5} simulation_${SIM_NAME5}.F90)
TARGET_LINK_LIBRARIES(
  ${EXECUTABLE5} 
  sll_simulations_sequential
  ${SIM_LIBS}
  hex_hermite
  sll_box_splines
  sll_hex_meshes)




# Installing simulations
INSTALL( TARGETS ${EXECUTABLE1} DESTINATION ${INSTALL_DIR} )
INSTALL( TARGETS ${EXECUTABLE2} DESTINATION ${INSTALL_DIR} )
INSTALL( TARGETS ${EXECUTABLE3} DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   ${INPUT_FILE}  README  DESTINATION ${INSTALL_DIR} )


# Ctest
# IF(BUILD_TESTING)
#   FILE( COPY ${INPUT_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#   ADD_TEST( NAME ${TEST_NAME1} COMMAND ${EXECUTABLE1} )
#   SET_TESTS_PROPERTIES(${TEST_NAME1} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
# ENDIF(BUILD_TESTING)

ENDIF (HDF5_ENABLED)
