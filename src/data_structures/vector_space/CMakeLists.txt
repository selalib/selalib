add_library( sll_vector_space STATIC
             sll_m_vector_space_base.F90
             sll_m_vector_space_real_array_1d.F90
             sll_m_vector_space_real_array_2d.F90
             sll_m_vector_space_real_array_3d.F90
             sll_m_vector_space_c1_block.F90
             sll_m_vector_space_c1_block_new.F90 )

target_link_libraries( sll_vector_space
                       sll_errors 
                       sll_working_precision)
# Tests
if( BUILD_TESTING )
  add_subdirectory( testing )
endif( BUILD_TESTING )
