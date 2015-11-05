SET(CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS_RELEASE})
ADD_LIBRARY(nufft2d STATIC dirft2d.f next235.f nufft2df90.f dfftpack.f)
ADD_EXECUTABLE(test_nufft2d nufft2d_demof90.f)
TARGET_LINK_LIBRARIES(test_nufft2d nufft2d)
