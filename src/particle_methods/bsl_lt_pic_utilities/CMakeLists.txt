
add_library( sll_bsl_lt_pic_utilities STATIC
             bsl_lt_pic_4d_utilities.F90 )

target_link_libraries( sll_bsl_lt_pic_utilities
  sll_pic_accumulators
  sll_assert
  sll_meshes
  sll_memory
  sll_timer
  sll_working_precision
  sll_file_io )

#  sll_bsl_lt_pic_particles

######## todo: OLD doc, please update

# ALH - <<test_lt_pic_4d_bsl_remap>> BSL remapping unit test
#
# source: [[file:unit_test_lt_pic_4d_bsl_remap.F90]]
# configuration, debug mode: [[elisp:(compile "cd ../../build && cmake -DCMAKE_BUILD_TYPE=Debug -DHDF5_PARALLEL_ENABLED=ON -DCMAKE_VERBOSE_MAKEFILE=TRUE -DOPENMP_ENABLED=OFF ../src")]]
# configuration, release mode: [[elisp:(compile "cd ../../build && cmake -DCMAKE_BUILD_TYPE=Release -DHDF5_PARALLEL_ENABLED=ON -DCMAKE_VERBOSE_MAKEFILE=FALSE -DOPENMP_ENABLED=OFF ../src")]]
# compile only: [[elisp:(compile "cd ../../build && make -j$NCPU test_lt_pic_4d_bsl_remap")]]
# single ltp run: [[elisp:(compile "cd ../../build && make -j $NCPU test_lt_pic_4d_bsl_remap && bin/test_lt_pic_4d_bsl_remap ltp")]]
# single ltp_bsl run: [[elisp:(compile "cd ../../build && make -j $NCPU test_lt_pic_4d_bsl_remap && bin/test_lt_pic_4d_bsl_remap ltp_bsl")]]
# debug: [[elisp:(gdb "gdb --annotate=3 ../../build/bin/test_lt_pic_4d_bsl_remap")]]
# clean: [[elisp:(compile "cd ../../build && make clean")]]
# ctest run: [[elisp:(compile "cd ../../build && ctest -R lt_pic_4d_bsl_remap")]]

#ADD_EXECUTABLE(test_lt_pic_4d_bsl_remap unit_test_lt_pic_4d_bsl_remap.F90)
#target_link_libraries(test_lt_pic_4d_bsl_remap
#  sll_pic_utilities
#  sll_pic_initializers
#  sll_pic_particles
#  sll_memory
#  sll_timer
#  sll_constants
#  sll_working_precision
#  sll_assert
#  )

#ADD_TEST(NAME test_lt_pic_4d_bsl_remap COMMAND test_lt_pic_4d_bsl_remap ltp_bsl)
#SET_TESTS_PROPERTIES(test_lt_pic_4d_bsl_remap PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20)

