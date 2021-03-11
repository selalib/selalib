add_library(sll_pic_weights STATIC sll_m_control_variate.F90)

target_link_libraries(sll_pic_weights sll_memory sll_working_precision
                      sll_particle_groups sll_distribution_function)

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
