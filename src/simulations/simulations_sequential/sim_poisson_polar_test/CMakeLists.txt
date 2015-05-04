set(  EXECUTABLE unit_test_polar                   )
set(   SRC_FILES unit_test_polar.F90               )
set(    DEP_LIBS sll_poisson_solvers               )
set( INSTALL_DIR ${INSTALL_DIR}/poisson_polar_test )

IF(MPI_ENABLED)
ADD_EXECUTABLE( ${EXECUTABLE} ${SRC_FILES} )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} ${DEP_LIBS} )

INSTALL( TARGETS ${EXECUTABLE} DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   README        DESTINATION ${INSTALL_DIR} )
ENDIF(MPI_ENABLED)

