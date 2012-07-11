add_library(sll_weno STATIC sll_WENO_node_interpolation.F90)

set(DEPENDENCIES sll_memory sll_assertion sll_working_precision)
add_library(WENO_interp STATIC WENO_interp.F90)
add_dependencies(WENO_interp ${DEPENDENCIES})
add_library(WENO_recon STATIC WENO_recon.F90)
add_dependencies(WENO_recon ${DEPENDENCIES})

add_executable(test_WENO unit_test.F90)
target_link_libraries(test_WENO
        WENO_interp
	WENO_recon
	sll_utilities
	sll_constants
	sll_assertion
	sll_memory
)
