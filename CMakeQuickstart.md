/I wrote this thanks to Trilinos CMake Quickstart/

Pierre

------------------------------------------------------------------------------
                         SeLaLib CMake Quickstart
------------------------------------------------------------------------------

A) Getting set up to use CMake
------------------------------

(*) Installing a binary release of CMake [Recommended for casual SeLaLib
    users]:

  Download and install the binary (currently version 2.8 is required) from:

    http://www.cmake.org/cmake/resources/software.html


(*) Installing CMake from source [Recommended for SeLaLib developers and
    experienced SeLaLib users]:

   wget http://www.cmake.org/files/v2.8/cmake-2.8.10.2.tar.gz
   tar zxvf cmake-2.8.10.2.tar.gz

   cd cmake-2.8.10.1; {
       ./bootstrap
       make
       make install
   }; cd -

  This will result in cmake and related CMake tools being installed in
  /usr/local/bin.

  NOTE: you probably must use sudo to install in the default /usr/local/bin.

B) Getting Help
---------------

(*) Finding CMake help at the website:

    http://www.cmake.org

(*) Building CMake help locally:

  $ cmake --help-full cmake.help.html

  (Open your web browser to the file cmake.help.html)

C) Configuring (Makefile Generator)
-----------------------------------

(*) Setting up a build directory:

    $ mkdir build
    $ cd build

  NOTE: You can create a build directory from any location you would like.  It
  can be a sub-directory of the SeLaLib base source directory or anywhere
  else.

  NOTE: If you mistakenly try to configure for an in-source build
  (e.g. with 'cmake .') you will get an error message and instructions
  on how to resolve the problem by deleting the generated
  CMakeCache.txt file (and other generated files) and then directions
  on how to create a different build directory as shown above.

(*) Basic configuration of SeLaLib:

  a) [Recommended] Create a 'do-configure' script such as:

        EXTRA_ARGS=$@
        
        cmake \
          -D CMAKE_BUILD_TYPE:STRING=Debug \
          $EXTRA_ARGS \
          ${SLL_HOME}

      and then run it with:

        ./do-configure [OTHER OPTIONS] -D<PACKAGE>_ENABLED=ON

      where <PACKAGE> FFTW, etc. and SLL_HOME is set to the
      selalib source base directory (or your can just give it explicitly).
      (ex: selalib/src)

  b) [Recommended] Create a CMake file fragment and point to it.

    Create a do-configure script like:

        EXTRA_ARGS=$@
        
        cmake \
          -D OPTIONS_FILE:FILEPATH=MyConfigureOptions.cmake \
          $EXTRA_ARGS \
          ${SLL_HOME}
       
    where MyConfigureOptions.cmake might look like:

      SET(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
      SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "" FORCE)
      SET(MPI_Fortran_COMPILER "mpiifort" CACHE BOOL "" FORCE) #Hydra build
      ...

    Using a configuration fragment file allows for better reuse of configure
    options across different configure scripts and better version control of
    configure options.

    NOTE: You can actually pass in a list of configuration fragment files
    which will be read in the order they are given.

    NOTE: If you do not use 'FORCE' shown above, then the option can be
    overridden on the cmake command line with -D options.  Also, if you don't
    use 'FORCE' then the option will not be set if it is already set in the
    case (e.g. by another configuration fragment file prior in the list).

  c) Using ccmake to configure:

    $ ccmake $SLL_HOME

  d) Using the QT CMake configuration GUI:

    On systems where it is installed, the CMake GUI
    can be a nice way to configure SeLaLib if you are a user.  To make your
    configuration easily repeatable, you might want to create a fragement file
    and just load it by setting OPTIONS_FILE (see above) in
    the GUI.

(*) Selecting the list of packages to enable:

  a) Configuring a package(s) along with all of the packages it can use:
  
      $ ./do-configure \
         -D HDF5_ENABLED:BOOL=ON \
         -D HDF5_PARALLEL_ENABLED:BOOL=ON \
         -D MPI_ENABLED:BOOL=ON \
         -D BUILD_TESTS:BOOL=ON
  
  b) Configuring SeLaLib to test all effects of changing a given package(s):
  
    $ ./do-configure \
       -D HDF5_ENABLE:BOOL=ON \
       -D FFTW_ENABLED:BOOL=OFF
  

(*) Selecting compiler and linker options:

  NOTE: The SeLaLib CMake build system will set up default compile options
  for GCC ('GNU') in development mode on order to help produce portable code.

  a) Configuring to build with default debug or release compiler flags:
  
    To build a debug version, pass into 'cmake':
  
      -D CMAKE_BUILD_TYPE:STRING=Debug
  
    This will result in default debug flags getting passed to the compiler.
  
    To build a release (optimized) version, pass into 'cmake':
  
      -D CMAKE_BUILD_TYPE:STRING=Release
  
    This will result in optimized flags getting passed to the compiler.
  
  b) Adding arbitrary compiler flags but keeping other default flags:

    To append arbitrary compiler flags that apply to all build types,
    configure with:

      -DCMAKE_Fortran_FLAGS:STRING="<EXTRA_COMPILER_OPTIONS>"

    where <EXTRA_COMPILER_OPTIONS> are your extra
    compiler options like "-DSOME_MACRO_TO_DEFINE -funroll-loops -fopenmp".  These
    options will get appened to other internally defined compiler option and
    therefore override them.

    NOTES:
 
    1) Setting CMAKE_Fortran_FLAGS with override but will not replace any other
    internally set flags in CMAKE_Fortran_FLAGS defined by the SeLaLib CMake
    system.  To get rid of these default flags, see below.

    2) For Fortran compiler, CMake passes compiler options to the compiler in the order:

      CMAKE_Fortran_FLAGS   CMAKE_Fortran_FLAGS_<CMAKE_BUILD_TYPE>

    where <CMAKE_BUILD_TYPE> = Debug or Release.  THEREFORE: The options in 
    CMAKE_Fortran_FLAGS_<CMAKE_BUILD_TYPE> come after and override those in 
    CMAKE_Fortran_FLAGS!.

    3) CMake defines default CMAKE_Fortran_FLAGS_<CMAKE_BUILD_TYPE> values that
    are overridden by the SeLaLib CMake build system for GCC ("GNU")
    compilers.
    This is mostly to provide greater control over the SeLaLib development
    environment.  This means that users setting the CMAKE_Fortran_FLAGS will
    *not* override the internally set debug or release flags in
    CMAKE_Fortran_FLAGS_<CMAKE_BUILD_TYPE> which come after on the compile
    line.  Therefore, setting CMAKE_Fortran_FLAGS should only be used for
    options that will not get overridden by the internally-set debug or
    release compiler flags in CMAKE_Fortran_FLAGS_<CMAKE_BUILD_TYPE>.  However,
    setting CMAKE_Fortran_FLAGS will work well for adding extra compiler
    defines (e.g. -DSOMETHING) for example.

    WARNING: Any options that you set through the cache varible
    CMAKE_Fortran_FLAGS_<CMAKE_BUILD_TYPE> (where <CMAKE_BUILD_TYPE> = Debug or
    Release) will get overridden in the SeLaLib CMake system for GNU
    compilers in development mode so don't try to manually set
    CMAKE_Fortran_FLAGS_<CMAKE_BUILD_TYPE>!
  
  c) Appending link flags to every executable:
  
    To pass in extra linker flags that are not
    libraries, use the built-in CMake variable CMAKE_EXE_LINKER_FLAGS.
  
(*) Configuring SeLaLib for MPI support:

  MPI is enabled by default, do disable the support you must minimally:

    -D MPI_ENABLED:BOOL=OFF

  a) Configuring build using MPI compiler wrappers:

    The MPI compiler wrappers are turned on by default.  There is built-in
    logic that will try to find the right compiler wrappers.  However, you can
    specifically select them by setting:

      -D MPI_[C,CXX_Fortran]_COMPILER:FILEPATH="exec_name"

        The name of the MPI C/C++/Fortran compiler wrapper executable.
        If this is just the name of the program it will be looked for
        in ${MPI_BIN_DIR} and in other standard locations with that name.
        If this is an absolute path, then this will be used as
        CMAKE_[C,CXX,Fortran]_COMPILER to compile and link code.

  b) Setting up to run MPI programs

    In order to use the ctest program to run MPI tests, you must set the mpi
    run command and the options it takes.  The built-in logic will try to find
    the right program and options but you will have to override them in many
    cases.

    MPI test and example executables are run as:

        ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
	               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXEC_NAME}
	               ${MPIEXEC_POSTFLAGS} ${ARGS})

    where ARGS, and PROCS are specific to the test being run.

    The test-independent MPI arguments are:

      -D MPI_EXEC:FILEPATH="exec_name"

        The name of the MPI run command (e.g. mpirun, mpiexec) that is used to
        run the MPI program.  This can be just the name of the program in which
        case the full path will be looked for in ${MPI_BIN_DIR} as described
        above.  If it is an absolute path, it will be used without question.

      -D MPIEXEC_MAX_NUMPROCS:STRING=4

        The maximum number of processes to allow when setting up and running
        MPI test and example executables.  The default is set to '4' and only
        needs to be changed when needed or desired.

      -D MPIEXEC_NUMPROC_FLAG:STRING=-np

        The command-line option just before the number of processes to use
        <NP>.  The default value is based on the name of ${MPI_EXEC}.

      -D MPIEXEC_POSTFLAGS:STRING="arg1 arg2 ... argn"

        Other command-line arguments that must come *after* the numprocs
        argument.  The default is empty "".


(*) Configuring SeLaLib for OpenMP support:

  To enable OpenMP support you must set

    -D OPENMP_ENABLED:BOOL=ON

  Note that if you enable OpenMP directly through a compiler option
  (e.g., -fopenmp), you will NOT enable OpenMP inside SeLaLib source code.

(*) Building shared libraries:

    -D BUILD_SHARED_LIBS:BOOL=ON

  NOTE: The above option will result in all shared libraries to be build on
  all systems (i.e. *.so on Unix/Linux systems, *.dylib on Mac OS X).


(*) Building static libraries and exectables:

   To build static libraries, turn off the shared library support:
  
    -D BUILD_SHARED_LIBS:BOOL=OFF


(*) Enabling support for optional BLAS/LAPACK:

<TODO>
  The headers, libraries, and library directories can then be specified with
  the input cache variables:

      Example:

        -D BLAS_LIBRARY_NAMES:STRING="blas;gfortran"

    BLAS_LIBRARY_DIRS:PATH: The list of directories where the
      library files can be found.

      Example:

        -D BLAS_LIBRARY_DIRS:PATH=/usr/local/blas
</TODO>


(*) Disabling tentatively enabled TPLs:

    -D BLAS_ENABLE_<TPLNAME>:BOOL=OFF

  NOTE: Some libraries in SeLaLib are always tentatively enabled (e.g. Lapack)
  and if all of the components are found (e.g. headers and libraries) then support will be enabled,
  otherwise it will be disabled.  
  It is possible that the enable process for the library may pass, but the
  library may not work correctly on the given platform.  In this case, one would
  also want to explicitly disable the TPL as shown above.


(*) Getting verbose output from the makefile:

    $ ./do_configure -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE


(*) Reconfiguring from scratch
  
      $ rm CMakeCache.txt ; find . -name CMakeFiles -exec rm -rf {} \;
      $ ./do-configure
  
    NOTE: Removing the CMakeCache.txt file is often needed when removing
    variables from the configure line.  Removing the CMakeFiles directories is
    needed if there are changes in some CMake modules or the CMake version
    itself.

(*) Viewing configure errors:

  Configure time errors are shown in the file:

      $BUILD_BASE_DIR/CMakeFiles/CMakeError.log

D) Building (Makefile generator)
--------------------------------

(*) Building all targets:

     $ make [-jN]

   (where N is the number of processes to use)


(*) Discovering what targets are available to build after configuration:

     $ make help


(*) Building with verbose output without reconfiguring:

    $ make [<SOME_TARGET>] VERBOSE=1


(*) Every executables are build in $BUILD_DIR/bin directory.


E) Testing with CTest
---------------------


(*) [Recommended] Testing using 'ctest'

    $ ctest -j4

  (see output in Testing/Temporary/LastTest.log)

  NOTE: The -jN argument allows CTest to use more processes to run
  tests but will not exceed the max number of processes specified at
  configure time.

  See detailed test output with:

    $ ctest -j4 -VV


(*) Only running tests for a single package

  Running a single test:

    $ ctest -j4 -R '<SLL_TEST>'

  (e.g. SLL_TEST = memory, collective, etc.)
  (see output in Testing/Temporary/LastTest.log)

(*) Running a single test with full output to the console:

    $ ctest -R ^FULL_TEST_NAME$ -VV

(*) Running memory checking:

  To run the memory tests for just a single package, from the *base* build
  directory, run:

    $ ctest -R '^<TEST>_' -T memcheck

  (see the detailed output in
  ./Testing/Temporary/LastDynamicAnalysis_DATE_TIME.log)

  NOTE: If you try to run memory tests from any subdirectories, that does not
  seem to work.  You have to run them from the base build directory and then
  use -R '^<TEST>_' with ctest in order to run your packages tests.


(*) Testing using 'make test'

    $ make test

  NOTE: This is equivalent to just running 'ctest' without any arguments.



F) Installing
-------------


(*) Setting the install prefix at configure time

    $ ./do-configure \
      -D CMAKE_INSTALL_PREFIX:PATH=$HOME/local \
      -D BUILD_PACKAGE:BOOL=ON

  NOTE: The script 'do-configure' is just a simple shell script that calls
  CMake as shown above.


(*) Installing after configuration

    $ make install

    (will build all of the targets needed before the install)


(*) Uninstall

    $ make uninstall



G) Packaging
------------


(*) Creating a tarball of the source tree:

   $ make package_source

   NOTE: The above command will tar up *everything* in the source tree (except
   for files explicitly excluded in the CMakeLists.txt files) so make sure
   that you start with a totally clean source tree before you do this.  Or,
   you could build Doxygen documentation first and then tar up SeLaLib and
   that would give you the source with Doxygen documentation.

   NOTE: You can control what gets put into the tarball by setting the cache
   variable CPACK_SOURCE_IGNORE_FILES when configuring with CMake.


H) Dashboard submissions
------------------------

You can use the extended CTest scripting system in SeLaLib to submit
package-by-package build, test, coverage, memcheck results to the dashboard.

First, configure as normal but add the build and test parallel levels with:

  $ ./do-configure -DCTEST_BUILD_FLAGS:STRING=-j4 -DCTEST_PARALLEL_LEVEL:STRING=4 \
    [OTHER OPTIONS]

Then, invoke the build, test and submit with:

  $ make Experimental

This invokes the advanced CTest script to do an experimental build
for all of the tests that you are enabled.  

*/