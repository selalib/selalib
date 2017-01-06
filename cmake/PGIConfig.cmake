IF(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  
  FIND_PROGRAM(PYTHON_EXECUTABLE NAMES python)
  FIND_PROGRAM(CPP_EXECUTABLE NAMES cpp)
  SET(PREPROCESS_SCRIPT "${CMAKE_SOURCE_DIR}/python/cut_long_lines.py")
  
  FUNCTION(CUT_LONG_LINES _FILE)
     MESSAGE(STATUS "Cut long lines of file ${_FILE}.F90")

     EXECUTE_PROCESS(COMMAND ${CPP_EXECUTABLE} "-Iinclude" "-E" "-w" "-P" 
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     INPUT_FILE "${CMAKE_SOURCE_DIR}/${_FILE}.F90"
                     OUTPUT_FILE "${CMAKE_SOURCE_DIR}/${_FILE}_pgi.F90")

     EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} ${PREPROCESS_SCRIPT} 
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

     SET(PGI TRUE CACHE BOOL "TRUE if using PGI fortran compiler")

  ENDFUNCTION(CUT_LONG_LINES)
  
  CUT_LONG_LINES(src/splines/sll_m_cubic_splines)
  CUT_LONG_LINES(src/quadrature/sll_m_gauss_legendre_integration)
# CUT_LONG_LINES(fft/sllfft/sll_fft)
# CUT_LONG_LINES(pic_utilities/sll_pic_utilities)
# CUT_LONG_LINES(pic_utilities/unit_test_particle_sort)
# CUT_LONG_LINES(pic_particle_initializers/sll_particle_init2D)
# CUT_LONG_LINES(pic_particle_initializers/sll_particle_init4D)
# CUT_LONG_LINES(interpolators/sll_cubic_spline_interpolator_1d)
# CUT_LONG_LINES(interpolators/sll_cubic_spline_interpolator_1d_nonuniform)
# CUT_LONG_LINES(interpolators/sll_cubic_spline_interpolator_2d)
# CUT_LONG_LINES(interpolators/sll_arbitrary_degree_spline_interpolator_1d)
# CUT_LONG_LINES(interpolators/sll_arbitrary_degree_spline_interpolator_2d)



ENDIF(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
