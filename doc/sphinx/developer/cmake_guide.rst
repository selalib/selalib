=====
CMake
=====
Add subdirectory
^^^^^^^^^^^^^^^^
To add a subdirectory to the build, modify the main CMakeList.txt in the *src* directory and add the following in the (*SUBDIRECTORY*) section ::

 add_subdirectory(<FOLDER_NAME>)

Example :

# src/CMakeList.txt ::

 add_subdirectory(fft)

Add module
^^^^^^^^^^
A module is a library. To add library to the current CMakeList.txt just add the lines ::

 add_library(<LIBRARY_NAME> STATIC <SOURCE.F90>)
 target_link_libraries(<LIBRARY_NAME> <LIST_OF_DEPENDENCIES>)

The <LIST_OF_DEPENDENCIES> is a collection of library that is required to link with build sll_fft.

Example :

# src/fft/CMakeList.txt::

 add_library(sll_fft STATIC sll_fft.F90)
 target_link_libraries(sll_fft sll_working_precision sll_memory sll_assert)

Add program
^^^^^^^^^^^

A program is an executable. In SeLaLib the executable are in build/bin/.
To add executable to the cmake tree add this lines in the current CMakeList::

 add_executable(<EXECUTABLE_NAME> <FILE_NAME.F90>)
 target_link_libraries(<EXECUTABLE_NAME> <LIBRARY_NAME>)

Example :

# src/fft/CMakeList.txt::

 add_executable(test_fft unit_test.F90)
 target_link_libraries(test_fft sll_fft) 

Add test
^^^^^^^^

A test is an executable define by ``add_executable`` that is automaticaly run by ``make test``.
Add the following lines to CMakeList.txt::

 add_test(NAME <TEST_NAME> COMMAND <EXECUTABLE_NAME>)
 set_tests_properties( <TEST_NAME> PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20)

where ``<TEST_NAME>`` is the name that appears when you run ``make
test`` and ``<EXECUTABLE_NAME>`` is the name of the executable that
you defined. The string "PASSED" could be changed and if needed,
increase the timeout value (seconds).

Example :

# src/fft/CMakeList.txt ::

  add_executable(test_fft unit_test.F90)
  target_link_libraries(test_fft sll_fft) 
  add_test(NAME fft_unit_test COMMAND test_fft)
  set_tests_properties( fft_unit_test PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20)

# output ::

  cd /selalib/build
  make test
  Running tests...
  Test project /Users/samuel/selalib/build
       Start 1: fft_unit_test
  1/12 Test #1: fft_unit_test ..............................   Passed    0.13 sec
  ...

Run ``make test`` in the buid directory is equivalent to run ``ctest``

Running individual Tests
^^^^^^^^^^^^^^^^^^^^^^^^

To determine what tests are available, you can always run ::

 ctest -N

which will display the list of tests but not actually run them. ::

 Test project /Users/samuel/selalib/build
  Test  #1: memory
  Test  #2: assert
  Test  #3: constants
  Test  #4: utilities
  Test  #5: low_level_file_io
  Test  #6: timer
  Test  #7: tridiagonal
  Test  #8: newton_raphson
  Test  #9: splines
  Test #10: integration
  Test #11: fft
  Test #12: collective
  Test #13: remap
  Test #14: WENO
  Test #15: interpolators
  Test #16: mapped_meshes
  Test #17: fields
  Test #18: ode_solvers
  Test #19: distribution_function
  Test #20: advection_field
  Test #21: poisson_solvers

 Total Tests: 21
 

The way of specifying tests is using explicit test number option -I ::

 ctest -I 3,5

will run tests ::

  Test project /Users/samuel/selalib/build
      Start 3: constants
  1/3 Test #3: constants ........................   Passed    0.00 sec
      Start 4: utilities
  2/3 Test #4: utilities ........................   Passed    9.09 sec
      Start 5: low_level_file_io
  3/3 Test #5: low_level_file_io ................   Passed    0.12 sec
  
  100% tests passed, 0 tests failed out of 3
  
  Total Test time (real) =   9.38 sec

Run ``ctest -I 3,3`` to run only the test #3.

If we now run ::

 ctest -R field

We will only see tests that contain string field ::

  Test project /Users/samuel/selalib/build
  
      Start 17: fields
  1/2 Test #17: fields ...........................   Passed    0.00 sec
      Start 20: advection_field
  2/2 Test #20: advection_field ..................   Passed    0.01 sec
  
  100% tests passed, 0 tests failed out of 2
  
  Total Test time (real) =   0.15 sec

We can also omit tests using -E, for example ::

 ctest -E field

will produce ::

  Test project /Users/samuel/selalib/build
        Start  1: memory
   1/19 Test  #1: memory ...........................   Passed    0.81 sec
        Start  2: assert
   2/19 Test  #2: assert ...........................   Passed    0.09 sec
        Start  3: constants
   3/19 Test  #3: constants ........................   Passed    0.00 sec
        Start  4: utilities
   4/19 Test  #4: utilities ........................   Passed    9.00 sec
        Start  5: low_level_file_io
   5/19 Test  #5: low_level_file_io ................   Passed    0.44 sec
        Start  6: timer
   6/19 Test  #6: timer ............................   Passed    2.00 sec
        Start  7: tridiagonal
   7/19 Test  #7: tridiagonal ......................   Passed    0.02 sec
        Start  8: newton_raphson
   8/19 Test  #8: newton_raphson ...................   Passed    0.02 sec
        Start  9: splines
   9/19 Test  #9: splines ..........................   Passed    0.01 sec
        Start 10: integration
  10/19 Test #10: integration ......................   Passed    0.00 sec
        Start 11: fft
  11/19 Test #11: fft ..............................   Passed    0.13 sec
        Start 12: collective
  12/19 Test #12: collective .......................   Passed    1.35 sec
        Start 13: remap
  13/19 Test #13: remap ............................   Passed    1.83 sec
        Start 14: WENO
  14/19 Test #14: WENO .............................   Passed    0.00 sec
        Start 15: interpolators
  15/19 Test #15: interpolators ....................   Passed    0.00 sec
        Start 16: mapped_meshes
  16/19 Test #16: mapped_meshes ....................   Passed    0.01 sec
        Start 17: ode_solvers
  17/19 Test #17: ode_solvers ......................   Passed    0.00 sec
        Start 18: distribution_function
  18/19 Test #18: distribution_function ............   Passed    0.02 sec
        Start 19: poisson_solvers
  19/19 Test #19: poisson_solvers ..................   Passed    0.00 sec
  
  100% tests passed, 0 tests failed out of 19

  Total Test time (real) =  15.78 sec
