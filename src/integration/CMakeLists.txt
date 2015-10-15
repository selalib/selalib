ADD_LIBRARY( sll_integration STATIC gauss_legendre.F90
                                    gauss_lobatto.F90
                                    rectangle.F90
                                    trapz.F90
                                    gauss_triangle.F90
                                    fekete.F90)

TARGET_LINK_LIBRARIES( sll_integration sll_splines
                                       sll_meshes
                                       sll_memory
                                       sll_utilities)
