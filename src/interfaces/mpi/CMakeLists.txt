if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions( -DHOMEBREW )
endif()

add_library( sll_mpi OBJECT sll_mpi.F90 )
