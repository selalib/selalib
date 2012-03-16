set(MESH_SRC sll_mesh_types.F90 geometry_functions.F90)
add_library(sll_mesh_types STATIC ${MESH_SRC})
set(MESH_LIBS sll_memory sll_working_precision sll_assert sll_utilities
	sll_constants sll_splines sll_diagnostics)
# add first test
set(SRC_TEST1 unit_test.F90)
add_executable(test_mesh_types ${SRC_TEST1})
target_link_libraries(test_mesh_types sll_mesh_types ${MESH_LIBS})

