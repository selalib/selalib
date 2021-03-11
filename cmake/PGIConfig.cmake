if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)

  set(PGI
      TRUE
      CACHE BOOL "TRUE if using PGI fortran compiler")
  find_program(PYTHON_EXECUTABLE NAMES python)
  find_program(FPP_EXECUTABLE NAMES gfortran)
  set(PREPROCESS_SCRIPT "${CMAKE_SOURCE_DIR}/python/cut_long_lines.py")

  function(CUT_LONG_LINES _FILE)
    message(STATUS "Cut long lines of file ${_FILE}.F90")

    execute_process(
      COMMAND ${FPP_EXECUTABLE} "-Iinclude" "-E" "-P"
              "${CMAKE_SOURCE_DIR}/${_FILE}.F90"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      OUTPUT_FILE "${CMAKE_SOURCE_DIR}/${_FILE}_pgi.F90")

    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${PREPROCESS_SCRIPT}
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

  endfunction(CUT_LONG_LINES)

  cut_long_lines(src/splines/sll_m_cubic_splines)
  cut_long_lines(src/quadrature/sll_m_gauss_legendre_integration)
  cut_long_lines(src/add_ons/multipatch/sll_m_cartesian_meshes_multipatch)
  cut_long_lines(src/field_solvers/maxwell_solvers/sll_m_maxwell_2d_pstd)
  cut_long_lines(src/field_solvers/maxwell_solvers/sll_m_maxwell_3d_pstd)
  cut_long_lines(
    src/field_solvers/maxwell_solvers_parallel/sll_m_maxwell_2d_periodic_cartesian_par
  )
  cut_long_lines(
    src/particle_methods/pic_opt/pic_opt_particle_initializers/sll_m_particle_initializers_2d
  )
  cut_long_lines(
    src/particle_methods/pic_opt/pic_opt_particle_initializers/sll_m_particle_initializers_4d
  )
  cut_long_lines(
    src/particle_methods/pic_opt/pic_opt_utilities/sll_m_pic_utilities)
  cut_long_lines(
    src/particle_methods/pic_basic/particle_groups/sll_m_particle_group_2d2v_lbf
  )
  cut_long_lines(src/interfaces/fft/sll_m_fft_sllfft)

endif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
