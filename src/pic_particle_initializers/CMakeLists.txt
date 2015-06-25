#HACK: THIS SHOULD BE MADE BETTER AND AT A LOWER LEVEL. Problem: when hdf5 is
# present and is parallel, we need to use the MPI compiler wrapper and this
# should propagate to all its clients, even if these were envisioned as
# sequential modules. Or, maybe this should be handled at the highest possible
# level: once the parallel hdf5 is found, all the library should be compiled 
# with the mpi wrapper compiler.
IF(HDF5_ENABLED AND HDF5_IS_PARALLEL)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_ENABLED AND HDF5_IS_PARALLEL)

# MLT_PIC conditional compilation - cf [[file:../pic_particle_types/CMakeLists.txt::MLT_PIC_OPTION]]

if(MLT_PIC)
  set(MLT_PIC_SRC_LIST
    mlt_pic_2d_init.F90 # [[file:mlt_pic_2d_init.F90]]
    )

  # ALH - <<test_mlt_pic_2d_init>>
  #
  # compile: [[elisp:(compile "cd ../../../build && make test_mlt_pic_2d_init && ctest -R mlt_pic_2d_init")]]

  ADD_EXECUTABLE(test_mlt_pic_2d_init unit_test_mlt_pic_2d_init.F90) # [[file:unit_test_mlt_pic_2d_init.F90]]
  TARGET_LINK_LIBRARIES(test_mlt_pic_2d_init
    sll_pic_particles
    sll_assert
    sll_constants
    sll_logical_meshes
    sll_pic_initializers
    )
  ADD_TEST(NAME mlt_pic_2d_init COMMAND test_mlt_pic_2d_init)
  SET_TESTS_PROPERTIES(mlt_pic_2d_init PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20)
endif()

ADD_LIBRARY( sll_pic_initializers STATIC 
  sll_particle_init4D.F90 
  sll_particle_init2D.F90
  lt_pic_4d_init.F90 
  ${MLT_PIC_SRC_LIST}
  )

TARGET_LINK_LIBRARIES( sll_pic_initializers
  sll_pic_utilities
  sll_pic_particles
#  sll_visu_pic
  sll_meshes
  sll_working_precision
  sll_assert
  sll_memory
  sll_random_deviate_generators )

# [[file:unit_test_lt_pic_4d_init.F90]]
IF(BUILD_TESTING)
  ADD_EXECUTABLE(lt_pic_4d_init_tester unit_test_lt_pic_4d_init.F90  )

  TARGET_LINK_LIBRARIES(lt_pic_4d_init_tester
    sll_pic_initializers
    sll_pic_utilities
    #  sll_lt_pic_4d_utilities
    sll_timer
    )

  ADD_TEST(NAME lt_pic_4d_init COMMAND lt_pic_4d_init_tester)
  SET_TESTS_PROPERTIES(lt_pic_4d_init PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20)
ENDIF(BUILD_TESTING)
  
# Ctest
IF(MPI_ENABLED AND BUILD_TESTING)
  ADD_EXECUTABLE( test_pic_initializers unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_pic_initializers
    sll_pic_initializers
    sll_collective )
  ADD_TEST( NAME pic_initializers COMMAND test_pic_initializers )
  SET_TESTS_PROPERTIES( pic_initializers PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(MPI_ENABLED AND BUILD_TESTING)
