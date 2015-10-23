set(    SIM_NAME eul_vp_2d2v_cart_fv          )
set(   TEST_NAME "sim_${SIM_NAME}"            )
set(  EXECUTABLE "${TEST_NAME}"               )
set( INSTALL_DIR "${INSTALL_DIR}/${SIM_NAME}" )

file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )

ADD_EXECUTABLE( ${EXECUTABLE} sll_m_sim_${SIM_NAME}.F90 ${TEST_NAME}.F90 )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} 
  sll_sim_parallel
  sll_parallel_array_utilities
  sll_coordinate_transformations
  sll_poisson_parallel
  sll_point_to_point_communications
  sll_meshes
  sll_fields
  sll_mesh_calculus
  sll_quadrature
  sll_file_io_parallel
  sll_timer
  sll_collective)

INSTALL( TARGETS ${EXECUTABLE}                     DESTINATION ${INSTALL_DIR} )
INSTALL( FILES  ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

