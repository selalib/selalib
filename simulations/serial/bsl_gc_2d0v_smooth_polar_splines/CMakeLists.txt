set(DIR_NAME bsl_gc_2d0v_smooth_polar_splines)
set(INSTALL_DIR ${INSTALL_DIR}/${DIR_NAME})
set(SIM_NAME sim_${DIR_NAME})

file(GLOB INPUT_FILES "*.nml")
file(GLOB POST_PROC "*.py")

# -------------------------------------------------------------------------------

add_library(
  ${SIM_NAME}_obj OBJECT
  sll_m_diagnostics.F90 sll_m_point_charge.F90 sll_m_simulation_state.F90
  sll_m_time_integrator_base.F90 sll_m_time_integrator_implicit.F90
  sll_m_time_integrator_explicit.F90)

add_dependencies(${SIM_NAME}_obj ${SLL_LIB})
add_executable(${SIM_NAME} ${SIM_NAME}.F90 $<TARGET_OBJECTS:${SIM_NAME}_obj>)
target_link_libraries(${SIM_NAME} ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS})

add_executable(${SIM_NAME}_equilibrium ${SIM_NAME}_equilibrium.F90)
target_link_libraries(${SIM_NAME}_equilibrium ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS})

# -------------------------------------------------------------------------------

install(TARGETS ${SIM_NAME} DESTINATION ${INSTALL_DIR})
install(TARGETS ${SIM_NAME}_equilibrium DESTINATION ${INSTALL_DIR})
install(FILES ${INPUT_FILES} ${POST_PROC} DESTINATION ${INSTALL_DIR})
