set( DIR_NAME bsl_gc_2d0v_smooth_polar_splines )
set( INSTALL_DIR ${INSTALL_DIR}/${DIR_NAME}    )
set( SIM_NAME sim_${DIR_NAME} )

#-------------------------------------------------------------------------------

add_library( ${SIM_NAME}_obj OBJECT )

add_dependencies( ${SIM_NAME}_obj ${SLL_LIB} )

add_executable( ${SIM_NAME} ${SIM_NAME}.F90 $<TARGET_OBJECTS:${SIM_NAME}_obj> )

target_link_libraries( ${SIM_NAME} ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS} )

#-------------------------------------------------------------------------------

install( TARGETS ${SIM_NAME} DESTINATION ${INSTALL_DIR} )
