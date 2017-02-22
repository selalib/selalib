SET( PIC_SAMPLING_SRCS
  sll_m_particle_sampling.F90 )

SET( PIC_SAMPLING_DEPS
  sll_pic_weights
  sll_memory
  sll_errors
  sll_working_precision
  sll_particle_groups
  sobol
  prob
  sll_distribution_function )

if (NOT (FFT_LIB STREQUAL "SLLFFT"))
  SET ( PIC_SAMPLING_SRCS
    ${PIC_SAMPLING_SRCS}
    sll_m_particle_sampling_interface.F90)
endif()

ADD_LIBRARY( sll_pic_sampling STATIC 
  ${PIC_SAMPLING_SRCS})

TARGET_LINK_LIBRARIES( sll_pic_sampling 
  ${PIC_SAMPLING_DEPS})


# Tests
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
