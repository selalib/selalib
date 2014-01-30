
ADD_EXECUTABLE(landau EXCLUDE_FROM_ALL landau.F90)
TARGET_LINK_LIBRARIES(landau selalib dfftpack)

ADD_EXECUTABLE(polar_advection EXCLUDE_FROM_ALL polar_advection.F90)
TARGET_LINK_LIBRARIES(polar_advection selalib)

ADD_CUSTOM_TARGET(examples DEPENDS landau polar_advection)
