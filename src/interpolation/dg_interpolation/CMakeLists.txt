set(DG_INTERPOLATION_SRC_FILES sll_m_dg_interpolator_1d.F90)

set(DG_INTERPOLATION_DEPS sll_quadrature sll_constants sll_interpolators)

add_library(sll_dg_interpolation STATIC ${DG_INTERPOLATION_SRC_FILES})

target_link_libraries(sll_dg_interpolation ${DG_INTERPOLATION_DEPS})

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
