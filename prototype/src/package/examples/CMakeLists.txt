
ADD_EXECUTABLE(landau EXCLUDE_FROM_ALL landau.F90)
TARGET_LINK_LIBRARIES(landau selalib dfftpack)

ADD_EXECUTABLE(polar_advection EXCLUDE_FROM_ALL polar_advection.F90)
TARGET_LINK_LIBRARIES(polar_advection selalib)

SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

ADD_EXECUTABLE(parallel_advection EXCLUDE_FROM_ALL parallel_advection.F90)
TARGET_LINK_LIBRARIES(parallel_advection selalib)

ADD_EXECUTABLE(landau_parallel EXCLUDE_FROM_ALL landau_parallel.F90)
TARGET_LINK_LIBRARIES(landau_parallel selalib dfftpack)

ADD_EXECUTABLE(vlasov_poisson_parallel EXCLUDE_FROM_ALL vlasov_poisson_parallel.F90)
TARGET_LINK_LIBRARIES(vlasov_poisson_parallel selalib dfftpack)

ADD_EXECUTABLE(poisson_2d_parallel EXCLUDE_FROM_ALL poisson_2d_parallel.F90)
TARGET_LINK_LIBRARIES(poisson_2d_parallel selalib dfftpack)

ADD_CUSTOM_TARGET(examples DEPENDS landau 
                                   polar_advection 
                                   parallel_advection
                                   vlasov_poisson_parallel
                                   poisson_2d_parallel
                                   landau_parallel)
