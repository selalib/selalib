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
  CUT_LONG_LINES(src/add_ons/multipatch/sll_m_cartesian_meshes_multipatch)
  CUT_LONG_LINES(src/field_solvers/maxwell_solvers/sll_m_maxwell_2d_pstd)
  CUT_LONG_LINES(src/field_solvers/maxwell_solvers/sll_m_maxwell_3d_pstd)
  CUT_LONG_LINES(src/field_solvers/maxwell_solvers_parallel/sll_m_maxwell_2d_periodic_cartesian_par)
  CUT_LONG_LINES(src/particle_methods/pic_opt2d/pic_opt2d_particle_initializers/sll_m_particle_initializers_2d)
  CUT_LONG_LINES(src/particle_methods/pic_opt2d/pic_opt2d_particle_initializers/sll_m_particle_initializers_4d)
  CUT_LONG_LINES(src/particle_methods/pic_opt2d/pic_opt2d_utilities/sll_m_pic_utilities)
  CUT_LONG_LINES(src/particle_methods/pic_basic/particle_groups/sll_m_particle_group_2d2v_lbf)

ENDIF(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
