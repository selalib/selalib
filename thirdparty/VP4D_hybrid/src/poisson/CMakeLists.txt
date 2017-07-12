##########################################################################
# poisson directory
##########################################################################
set(files 
	jrk_testcase.F90
	jrk_model_def.F90 
	jrk_model.F90
	jrk_diagnostics_computation.F90
	jrk_element_contribution.F90)

set(mods 
	${CMAKE_BINARY_DIR}/modules/jrk_testcase.mod
	${CMAKE_BINARY_DIR}/modules/jrk_model_def.mod
	${CMAKE_BINARY_DIR}/modules/jrk_model.mod
	${CMAKE_BINARY_DIR}/modules/jrk_diagnostics_computation.mod
	${CMAKE_BINARY_DIR}/modules/jrk_element_contribution.mod)

add_library(jrk_poisson ${files})
target_link_libraries(jrk_poisson ${JOREK_LIBRARIES} ${SPM_LIBRARIES} ${MPI_LIBRARY})

install(TARGETS jrk_poisson DESTINATION lib)
install(FILES ${mods} DESTINATION include)
##########################################################################
