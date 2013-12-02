ADD_LIBRARY(sll_csl STATIC sll_csl.F90)

TARGET_LINK_LIBRARIES(sll_csl sll_distribution_function 
                              sll_interpolators 
                              ode_solvers
                              sll_scalar_field 
                              sll_file_io 
                              sll_utilities )

ADD_EXECUTABLE(test_csl unit_test.F90)
TARGET_LINK_LIBRARIES(test_csl sll_csl )

ADD_LIBRARY(sll_bgk_mesh STATIC sll_bgk_mesh.F90)
TARGET_LINK_LIBRARIES(sll_bgk_mesh sll_csl)

ADD_EXECUTABLE(vp_csl_bloc_unif vp_csl_bloc_unif.F90)
TARGET_LINK_LIBRARIES(vp_csl_bloc_unif sll_bgk_mesh)

ADD_EXECUTABLE(vp_classic_csl bgk_classic_csl.F90 
                              sll_classic_csl.F90 
                              contrib_rho.F90 )

TARGET_LINK_LIBRARIES(vp_classic_csl sll_bgk_mesh)

ADD_EXECUTABLE(vp_non_unif_mesh vp_non_unif_mesh.F90 
                                sll_classic_csl.F90 
                                contrib_rho.F90 )

TARGET_LINK_LIBRARIES(vp_non_unif_mesh sll_bgk_mesh)

#ADD_EXECUTABLE(sa_csl sa_csl.F90)
#TARGET_LINK_LIBRARIES(sa_csl sll_csl )