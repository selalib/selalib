IF(MPI_ENABLED)
# Executables must be linked against following Selalib libraries
set( DEP_LIBS sll_interpolators
              sll_fft
              sll_distribution_function
              sll_constants
              sll_file_io
  	          sll_utilities
              sll_poisson_solvers
              dfftpack )

# Installation directory
set( INSTALL_DIR "${INSTALL_DIR}/VP1D_various_programs" )

# NEW_PROGRAM macro: compiles, links and installs executable
MACRO( NEW_PROGRAM EXEC_NAME SOURCE_FILE )
  ADD_EXECUTABLE( ${EXEC_NAME} ${SOURCE_FILE} )
  TARGET_LINK_LIBRARIES( ${EXEC_NAME} ${SIM_LIBS} ${DEP_LIBS} )
  INSTALL( TARGETS ${EXEC_NAME} DESTINATION ${INSTALL_DIR} )
ENDMACRO( NEW_PROGRAM )

# New executables 
NEW_PROGRAM( VP1D_cart VP_uniform_cartesian_BSL.F90 )
IF( HDF5_ENABLED )
  NEW_PROGRAM( VP1D_deltaf_cart VP1D_deltaf_BSL.F90 )
  NEW_PROGRAM( VP1D_mm VP1D_mm.F90 )
ENDIF( HDF5_ENABLED )

# Copy additional files to installation directory
file( GLOB INPUT_FILES "*.nml"         )
file( GLOB   POST_PROC "*.dat" "*.gnu" )
INSTALL( FILES ${INPUT_FILES} README ${POST_PROC} DESTINATION ${INSTALL_DIR} )

ENDIF(MPI_ENABLED)
