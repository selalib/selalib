SET( INSTALL_DIR "${INSTALL_DIR}/post_processing_tools" )

IF(HDF5_ENABLED AND HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED AND MPI_MODULE_ENABLED)

   #===========================================================================
   # Add postprocessing code for KEEN
   IF(Fortran_COMPILER MATCHES IBM)
      SET_SOURCE_FILES_PROPERTIES( VPpostprocessing_drive_KEEN.f 
                                   PROPERTIES COMPILE_FLAGS -qfixed )
   ENDIF(Fortran_COMPILER MATCHES IBM)

   ADD_EXECUTABLE( VPpostprocessing_drive_KEEN VPpostprocessing_drive_KEEN.f )
   TARGET_LINK_LIBRARIES( VPpostprocessing_drive_KEEN dfftpack )
   INSTALL( TARGETS VPpostprocessing_drive_KEEN DESTINATION ${INSTALL_DIR} )
   #===========================================================================
   # Some misterious Python script for visualizing output data with matplotlib
   INSTALL( FILES vp_plots.py DESTINATION ${INSTALL_DIR} )
   #===========================================================================

ENDIF()
