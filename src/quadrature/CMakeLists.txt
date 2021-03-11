set(QUADRATURE_SRCS
    sll_m_gauss_lobatto_integration.F90 sll_m_rectangle_integration.F90
    sll_m_trapz_integration.F90 sll_m_fekete_integration.F90
    sll_m_primitives.F90)

if(PGI)
  list(APPEND QUADRATURE_SRCS sll_m_gauss_legendre_integration_pgi.F90)
else(PGI)
  list(APPEND QUADRATURE_SRCS sll_m_gauss_legendre_integration.F90)
endif(PGI)

add_library(sll_quadrature STATIC ${QUADRATURE_SRCS})

target_link_libraries(sll_quadrature sll_splines sll_meshes sll_memory
                      sll_utilities)
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
