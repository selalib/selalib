# DD version with dynamic halos
set(SIM_NAME_SLIM bsl_vp_3d3v_cart_dd_slim)
set(TEST_NAME_SLIM "sim_${SIM_NAME_SLIM}")
set(EXECUTABLE_SLIM "${TEST_NAME_SLIM}")

add_executable(
  ${EXECUTABLE_SLIM} sll_m_sim_${SIM_NAME_SLIM}.F90 ${TEST_NAME_SLIM}.F90 )

if(USE_FMEMPOOL)
  set(FMEMPOOL "fmempool")
  add_definitions(-DUSE_FMEMPOOL)
endif(USE_FMEMPOOL)

# build a C/C++ callable interface to simulation
enable_language(C)
enable_language(CXX)

include(FortranCInterface)

set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

add_library(
  sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface STATIC
  sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface.F90
  sll_m_sim_bsl_vp_3d3v_cart_dd_slim.F90 sll_m_sim_6d_utilities.F90)

target_link_libraries(sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface ${SLL_LIB}
                      ${SLL_EXT} ${SLL_DEPS})

install(TARGETS sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
        # DESTINATION ${INSTALL_DIR}
	)

# build a test that checks if the simulation's functions (exported via
# _interface module) can be called from a C++ program cf
# https://stackoverflow.com/questions/28344564/cmake-remove-a-compile-flag-for-a-single-translation-unit
macro(remove_flag_from_target _target _flag)
  get_target_property(_target_cxx_flags ${_target} COMPILE_OPTIONS)
  if(_target_cxx_flags)
    list(REMOVE_ITEM _target_cxx_flags ${_flag})
    set_target_properties(${_target} PROPERTIES COMPILE_OPTIONS
                                                "${_target_cxx_flags}")
  endif()
endmacro()

# cf.
# https://github.com/bast/python-cffi-demo/blob/exercise/cxx-fortran/CMakeLists.txt
# c executable
add_executable(test_cpp_interface test_cpp_interface.cpp)

set_target_properties(
  test_cpp_interface PROPERTIES LINKER_LANGUAGE Fortran CXX_STANDARD 17
                                CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
# In case of intel: remove compiler flags -nowarn -ipo-separate
remove_flag_from_target(test_cpp_interface -nowarn)
remove_flag_from_target(test_cpp_interface -ipo-separate)

# flag required for ifort target_link_options("-nofor_main")

# linked against fortran libraries
target_link_libraries(test_cpp_interface
	sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface ${MPI_CXX_LIBRARIES})
install(TARGETS test_cpp_interface DESTINATION ${INSTALL_DIR})

target_link_libraries(${EXECUTABLE_SLIM} sll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface
	${SLL_LIB} ${SLL_EXT} ${SLL_DEPS}
                      ${FMEMPOOL})
install(TARGETS ${EXECUTABLE_SLIM} DESTINATION ${INSTALL_DIR})

file(GLOB INPUT_FILES "*.nml")
file(GLOB REFERENCE_FILES "reffile_*.dat")
file(GLOB POST_PROC "*.dat" "*.gnu")

# TODO : add slim test, add slim overlap test
if(BUILD_TESTING)
  file(COPY ${INPUT_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY ${REFERENCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(ARGS ${CMAKE_CURRENT_SOURCE_DIR}/param_bsl_vp_3d3v_cart_dd_testing.nml)

  add_mpi_test(bsl_vp_3d3v_cart_dd_slim ${EXECUTABLE_SLIM} ${PROCS} ${ARGS})
  set_tests_properties(bsl_vp_3d3v_cart_dd_slim
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  add_mpi_test(test_cpp_interface test_cpp_interface ${PROCS} ${ARGS})
  set_tests_properties(test_cpp_interface PROPERTIES PASS_REGULAR_EXPRESSION
                                                     "works in cpp")
endif(BUILD_TESTING)
