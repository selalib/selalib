set(PIC_SAMPLING_SRCS sll_m_particle_sampling.F90)

set(PIC_SAMPLING_DEPS
    sll_pic_weights
    sll_memory
    sll_errors
    sll_working_precision
    sll_particle_groups
    sll_pic_weights
    sobol
    prob
    sll_distribution_function)

if(NOT (FFT_LIB STREQUAL "SLLFFT"))
  set(PIC_SAMPLING_SRCS ${PIC_SAMPLING_SRCS}
                        sll_m_particle_sampling_interface.F90)
endif()

add_library(sll_pic_sampling STATIC ${PIC_SAMPLING_SRCS})

target_link_libraries(sll_pic_sampling ${PIC_SAMPLING_DEPS})

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
