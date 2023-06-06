

# DD version with dynamic halos
set( SIM_NAME_SLIM bsl_vp_3d3v_cart_dd_slim )
set( TEST_NAME_SLIM "sim_${SIM_NAME_SLIM}" )
set( EXECUTABLE_SLIM "${TEST_NAME_SLIM}" )

ADD_EXECUTABLE( ${EXECUTABLE_SLIM}
                sll_m_sim_${SIM_NAME_SLIM}.F90
                ${TEST_NAME_SLIM}.F90
                sll_m_sim_6d_utilities.F90 )

IF(USE_FMEMPOOL)
  set( FMEMPOOL "fmempool" )
  add_definitions( -DUSE_FMEMPOOL )
ENDIF(USE_FMEMPOOL)

# build a C/C++ callable interface to simulation
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

include(FortranCInterface)

#TODO have this detected
SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
# # FortranCInterface_VERIFY(CXX)
# Message("fortran mpi libs " ${MPI_Fortran_LIBRARIES} )
# Message("fortran flags " ${MPI_Fortran_LINK_FLAGS} )
# Message("fortran cmake libs " ${CMAKE_Fortran_LIBRARIES} )
# # string(STRIP ${MPI_Fortran_LIBRARIES} strippedinput)
# string(REPLACE
#   "mpifh" "cxx" CXX_ADD_LIBRARIES_PRE
#   "${MPI_Fortran_LIBRARIES}")
# string(REPLACE
#   "fort" "cxx" CXX_ADD_LIBRARIES
#   "${CXX_ADD_LIBRARIES_PRE}")
# set ( CMAKE_CXX_LIBRARIES ${CMAKE_CXX_LIBRARIES} ${CXX_ADD_LIBRARIES})
# 
# # we need /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so here or similar
# Message("cxx libs " ${CMAKE_CXX_LIBRARIES} )

# our interface library
ADD_LIBRARY( sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
              STATIC
              sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface.F90
              sll_m_sim_bsl_vp_3d3v_cart_dd_slim.F90
              sll_m_sim_6d_utilities.F90
              )
TARGET_LINK_LIBRARIES( sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
  sll_working_precision
  sll_memory
  ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS}
  stdc++
  ${MPI_CXX_LIBRARIES}
  )

#INSTALL( TARGETS sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
#  # DESTINATION ${INSTALL_DIR_BASE}/sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface 
#)

# build a test that checks if the simulation's functions (exported via _interface module)
# can be called from a C++ program
# cf https://stackoverflow.com/questions/28344564/cmake-remove-a-compile-flag-for-a-single-translation-unit
macro(remove_flag_from_target _target _flag)
    get_target_property(_target_cxx_flags ${_target} COMPILE_OPTIONS)
    if(_target_cxx_flags)
        list(REMOVE_ITEM _target_cxx_flags ${_flag})
        set_target_properties(${_target} PROPERTIES COMPILE_OPTIONS "${_target_cxx_flags}")
    endif()
endmacro()

# cf. https://github.com/bast/python-cffi-demo/blob/exercise/cxx-fortran/CMakeLists.txt
# c executable
add_executable(
  test_cpp_interface
  test_cpp_interface.cpp
    )
SET_TARGET_PROPERTIES(test_cpp_interface
  PROPERTIES
  LINKER_LANGUAGE Fortran
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  )
#In case of intel: remove compiler flags
# -nowarn -ipo-separate
remove_flag_from_target(test_cpp_interface -nowarn)
remove_flag_from_target(test_cpp_interface -ipo-separate)

# flag required for ifort
# target_link_options("-nofor_main")

# linked against fortran libraries
target_link_libraries(
  test_cpp_interface
  sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
  ${MPI_CXX_LIBRARIES}
#  "-nofor_main"
)
install(TARGETS test_cpp_interface DESTINATION ${INSTALL_DIR})

TARGET_LINK_LIBRARIES( ${EXECUTABLE_SLIM} ${SLL_LIB} ${SLL_EXT} ${SLL_DEPS} ${FMEMPOOL} )
install(TARGETS ${EXECUTABLE_SLIM} DESTINATION ${INSTALL_DIR})


file( GLOB INPUT_FILES "*.nml" )
file( GLOB REFERENCE_FILES "reffile_*.dat" )
file( GLOB POST_PROC "*.dat" "*.gnu" )

#INSTALL( FILES ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR_BASE} )


# TODO : add slim test, add slim overlap test
IF(BUILD_TESTING)
  FILE( COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  FILE( COPY ${REFERENCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
  SET( ARGS ${CMAKE_CURRENT_SOURCE_DIR}/param_bsl_vp_3d3v_cart_dd_testing.nml )

  ADD_MPI_TEST( bsl_vp_3d3v_cart_dd_slim ${EXECUTABLE_SLIM} ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( bsl_vp_3d3v_cart_dd_slim PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

  ADD_MPI_TEST( test_cpp_interface test_cpp_interface ${PROCS} ${ARGS} )
  SET_TESTS_PROPERTIES( test_cpp_interface PROPERTIES PASS_REGULAR_EXPRESSION "works in cpp" )
ENDIF(BUILD_TESTING)
