add_library(sll_boundary_condition_descriptors STATIC
            sll_m_boundary_condition_descriptors.F90)

target_link_libraries(sll_boundary_condition_descriptors sll_working_precision)
