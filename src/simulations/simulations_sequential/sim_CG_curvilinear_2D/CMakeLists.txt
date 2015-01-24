IF(HDF5_ENABLED)
  # Old library: Poisson solver on "Colella" curvilinear mesh, uses Mudpack
  ADD_LIBRARY( sll_mudpack_colella STATIC
               sll_mudpack_base.F90 sll_mudpack_colella.F90 )
  TARGET_LINK_LIBRARIES( sll_mudpack_colella
                         sll_working_precision sll_memory sll_assert mudpack )
  
  # Executable: unit-test for library above
  ADD_EXECUTABLE( unit_test_mudpack_colella test_mudpack_colella.F90 )
  TARGET_LINK_LIBRARIES( unit_test_mudpack_colella
                         sll_mudpack_colella sll_constants )
  
  # Executable: simulation
  ADD_EXECUTABLE( CG_Curvilinear_2D 
                  CG_Curvilinear_2D.F90 
                  CG_Curvilinear_2D_advection.F90
                  CG_Curvilinear_2D_diagnostic.F90
                  CG_Curvilinear_2D_structure.F90
                  CG_Curvilinear_Function.F90 )
  
  TARGET_LINK_LIBRARIES( CG_Curvilinear_2D
                         sll_constants
                         sll_poisson_solvers
                         sll_splines
                         sll_interpolators
                         sll_mudpack_colella
                         sll_fft
                         sll_utilities )
  
  # Installation of executables and other files
  set( INPUT_FILES CG_Curvilinear_data.txt          )
  set( INSTALL_DIR ${INSTALL_DIR}/CG_Curvilinear_2D )
  
  INSTALL( TARGETS unit_test_mudpack_colella CG_Curvilinear_2D
           DESTINATION ${INSTALL_DIR} )
  INSTALL( FILES   ${INPUT_FILES} README CG_2D_Curvilinear-exe.py
           DESTINATION ${INSTALL_DIR} )

ENDIF(HDF5_ENABLED)
