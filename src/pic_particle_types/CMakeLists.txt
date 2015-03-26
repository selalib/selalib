ADD_LIBRARY( sll_pic_particles STATIC 
  particle_representation.F90
  particle_group_2d.F90 
  particle_group_4d.F90
  lt_pic_4d_particle.F90
  lt_pic_4d_group.F90
  ${MLT_PIC_SRC_LIST}
  )

TARGET_LINK_LIBRARIES( sll_pic_particles
  sll_memory
  sll_working_precision
  sll_assert
  sll_meshes )

# <<MLT_PIC_OPTION>> - ALH - mlt_pic developments are only conditionnally compiled for the moment to allow unfinished
# program files to be exchanged with MCP without breaking the regression test procedure.  Activating MLT_PIC: cmake
# -DMLT_PIC=ON.  Deactivating MLT_PIC (just to make sure that Jenkins is OK): cmake -DMLT_PIC=OFF

OPTION(MLT_PIC "Activating mlt_pic developments and tests" OFF)

if(MLT_PIC)
  set(MLT_PIC_SRC_LIST
    mlt_pic_2d_particle.F90 # [[file:mlt_pic_2d_particle.F90]]
    mlt_pic_2d_group.F90 # [[file:mlt_pic_2d_group.F90]]
    )

  # ALH - <<test_mlt_pic_2d_particles>>
  #
  # compile: [[elisp:(compile "cd ../../../build && make test_mlt_pic_2d_particles && ctest -R mlt_pic_2d_particles")]]

  ADD_EXECUTABLE(test_mlt_pic_2d_particles unit_test_mlt_pic_2d.F90) # [[file:unit_test_mlt_pic_2d.F90]]
  TARGET_LINK_LIBRARIES(test_mlt_pic_2d_particles
    sll_pic_particles
    sll_assert
    sll_constants
    sll_meshes
    )
  ADD_TEST(NAME mlt_pic_2d_particles COMMAND test_mlt_pic_2d_particles)
  SET_TESTS_PROPERTIES(mlt_pic_2d_particles PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20)
endif()

# Ctest
IF(BUILD_TESTING)
  ADD_EXECUTABLE( test_pic_particles unit_test.F90 )
  TARGET_LINK_LIBRARIES( test_pic_particles sll_pic_particles )
  ADD_TEST( NAME pic_particles COMMAND test_pic_particles )
  SET_TESTS_PROPERTIES( pic_particles PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF(BUILD_TESTING)
