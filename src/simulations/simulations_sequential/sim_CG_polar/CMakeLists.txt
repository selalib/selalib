set(  EXECUTABLE CG_polar                )
set(   SRC_FILES CG_polar.F90 sll_advection_polar.F90 sll_polar_operators.F90 )
set( INPUT_FILES CG_data.txt             )
set( INSTALL_DIR ${INSTALL_DIR}/CG_polar )
set(    DEP_LIBS sll_constants
                 sll_fft
                 sll_file_io
                 sll_poisson_solvers
                 sll_splines
                 sll_timer
                 sll_utilities )

IF(MPI_ENABLED)
ADD_EXECUTABLE( ${EXECUTABLE} ${SRC_FILES} )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} ${DEP_LIBS} )

INSTALL( TARGETS ${EXECUTABLE}         DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   ${INPUT_FILES} README DESTINATION ${INSTALL_DIR} )
ENDIF(MPI_ENABLED)

