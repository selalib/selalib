set(CHARACTERISTICS_SRC_FILES
    sll_m_characteristics_1d_explicit_euler.F90
    sll_m_characteristics_1d_explicit_euler_conservative.F90
    sll_m_characteristics_1d_trapezoid.F90
    sll_m_characteristics_1d_trapezoid_conservative.F90
    sll_m_characteristics_1d_base.F90
    sll_m_characteristics_2d_verlet.F90
    sll_m_characteristics_2d_explicit_euler.F90
    sll_m_characteristics_2d_explicit_euler_conservative.F90
    sll_m_characteristics_2d_base.F90
    sll_m_euler_2d_hex.F90)

add_library(sll_characteristics STATIC ${CHARACTERISTICS_SRC_FILES})
target_link_libraries(sll_characteristics sll_utilities sll_interpolators)

# Ctests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
