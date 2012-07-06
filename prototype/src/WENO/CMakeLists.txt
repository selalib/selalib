
#add_library(sll_interpolators_1d STATIC ../interpolators/sll_interpolators_1d_base.F90)
#add_library(sll_weno STATIC sll_WENO_node_interpolation.F90)

add_library(WENO_interp STATIC WENO_interp.F90)
add_library(WENO_recon STATIC WENO_recon.F90)

add_executable(test_WENO_recon unit_test_weno_recon.F90)
add_executable(test_WENO_interp unit_test_weno_interp.F90)

target_link_libraries(test_WENO_recon
	WENO_recon
	sll_utilities
	sll_constants
	sll_assertion
	sll_memory
)

target_link_libraries(test_WENO_interp
        WENO_interp
        sll_utilities
        sll_constants
        sll_assertion
        sll_memory
)

