set(  EXECUTABLE polar_SL                )
set(   SRC_FILES polar_SL.F90            )
set( INSTALL_DIR ${INSTALL_DIR}/polar_SL )
set(    DEP_LIBS sll_interpolators sll_constants sll_file_io sll_utilities )

ADD_EXECUTABLE( ${EXECUTABLE} ${SRC_FILES} )
TARGET_LINK_LIBRARIES( ${EXECUTABLE} ${DEP_LIBS} )

INSTALL( TARGETS ${EXECUTABLE} DESTINATION ${INSTALL_DIR} )
INSTALL( FILES   README        DESTINATION ${INSTALL_DIR} )

