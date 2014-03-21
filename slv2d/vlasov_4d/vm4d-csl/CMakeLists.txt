ADD_LIBRARY(slv2d STATIC 
              diagnostiques.F90
              geometry.F90
              maxwell2dfdtd.F90
              poisson2dpp_seq.F90
              fft_module.F90
              precision.F90
              splinenn.F90
              splinepp.F90
              splinepx.F90
              splinepy.F90
              transpose.F90
              vlasov1d.F90
              vlasov2d.F90 
              vlasov4d_plot.F90)


ADD_EXECUTABLE(vm4d_transpose vm4d_transpose.F90 
                              ../init_functions.F90
               ../../../prototype/src/poisson_solvers/sll_mudpack_base.F90
               ../../../prototype/src/poisson_solvers/sll_mudpack_cartesian.F90)

TARGET_LINK_LIBRARIES(vm4d_transpose slv2d ${LIBS} 
                                           dfftpack
                                           mudpack
                                           ${LAPACK_LIBRARIES})
