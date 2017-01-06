SET(QUADRATURE_SRCS sll_m_gauss_lobatto_integration.F90
                    sll_m_rectangle_integration.F90
                    sll_m_trapz_integration.F90
                    #sll_m_gauss_triangle_integration.F90
                    sll_m_fekete_integration.F90
                    sll_m_primitives.F90)

IF(PGI)
  LIST(APPEND QUADRATURE_SRCS sll_m_gauss_legendre_integration_pgi.F90)
ELSE(PGI)
  LIST(APPEND QUADRATURE_SRCS sll_m_gauss_legendre_integration.F90)
ENDIF(PGI)

ADD_LIBRARY( sll_quadrature STATIC ${QUADRATURE_SRCS})

TARGET_LINK_LIBRARIES( sll_quadrature sll_splines
                                      sll_meshes
                                      sll_memory
                                      sll_utilities)
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
