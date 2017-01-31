FIND_PROGRAM(FORCHECK_EXECUTABLE
  NAMES forchk
  PATHS $ENV{FCKDIR}/bin /usr/bin /bin /usr/local/bin
  DOC "Performs a full static analysis of Fortran programs.")

IF (FORCHECK_EXECUTABLE)
  SET (FORCHECK_FOUND "YES")
ENDIF (FORCHECK_EXECUTABLE)

MARK_AS_ADVANCED(
  FORCHECK_FOUND
  FORCHECK_EXECUTABLE
  )

IF(FORCHECK_FOUND)

   # The result of the forcheck analysis will go into this directory
   set(FORCHECK_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/forcheck)
   file(MAKE_DIRECTORY ${FORCHECK_DIRECTORY})

   # Some files that are required for the analysis are here
   set(FORCHECK_FLB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/forcheck)

   # Forcheck library files for the external libraries
   set(FORCHECK_EXTERNAL_FLBS $ENV{FCKDIR}/share/forcheck/MPI_3.flb ${FORCHECK_FLB_DIR}/hdf5-1_8_9.flb)
   if(NOT EXISTS "$ENV{FCKDIR}/share/forcheck/MPI_3.flb")
     message(WARNING "Forcheck: Can't find MPI_3.flb.\n Most probably the Forcheck module is not loaded.\n Try to load it and rerun cmake.")
   endif()

   #ADD_CUSTOM_TARGET(forcheck "${FORCHECK_EXECUTABLE} -define DEBUG,GFORTRAN -I ${INCS} -l mylistfile -ff ${SRCS} $(FCKDIR)/share/forcheck/MPI.flb"  COMMENT "Forcheck the source code" VERBATIM)

ENDIF(FORCHECK_FOUND)


# adds a custom target for running the Forcheck analysis
# call this function at the end of the CMakeLists.txt
function(add_forcheck_target)
  # retriev the lists that were creatde by add_library:
  get_property(_fck_sources GLOBAL PROPERTY CPP_SOURCES)
  get_property(_fck_preproc_sources GLOBAL PROPERTY CPP_PREPROC_SOURCES)
  get_property(_fck_includes GLOBAL PROPERTY CPP_INCLUDES)

  # set up include flags for preprocessing
  set(incflags)
  foreach(i ${_fck_includes})
    set(incflags ${incflags} -I${i})
  endforeach()

  # list of compiler definitions to be used for preprocessing 
  foreach(_source ${_fck_sources})
     get_source_file_property(_defs "${_source}" COMPILE_DEFINITIONS)
     list(APPEND _fck_defines "${_defs}")
     get_filename_component(_dir ${_source} PATH)
     get_property( _defs  DIRECTORY ${_dir} PROPERTY COMPILE_DEFINITIONS)
     if(_defs)
       list(APPEND _fck_defines "${_defs}")
     endif()
  endforeach()
  if(_fck_defines)
     list(REMOVE_DUPLICATES _fck_defines)
  endif()
  set(_defflags)
  foreach(d ${_fck_defines})
    set(_defflags ${_defflags} -D${d})
  endforeach() 
 
  # Create custom commands for preprocessing the Fortran files
#   while(_fck_sources)
#      list(GET _fck_sources 0 _src)
#      list(GET _fck_preproc_sources 0 _preproc_src)
#      list(REMOVE_AT _fck_sources 0)
#      list(REMOVE_AT _fck_preproc_sources 0)
#      add_custom_command(OUTPUT "${_preproc_src}"
#          COMMAND ${CMAKE_Fortran_COMPILER} ${incflags} ${_defflags} ${preprocessor_only_flags} ${_src} | sed -f ${FORCHECK_FLB_DIR}/sedscript -f ${FORCHECK_FLB_DIR}/sedscript2  > ${_preproc_src}
#          #IMPLICIT_DEPENDS Fortran "${_source}"
#          DEPENDS "${_src}"
#          COMMENT "Preprocessing ${_src}"
#          VERBATIM
#        ) 
#      # the preprocessed file is piped through a sed script, to break up the long lines
#      #set_source_files_properties(${_preproc_src} PROPERTIES GENERATED TRUE)
#   endwhile()

  # group all preprocessing commands into one target
  #get_property(_fck_preproc_sources GLOBAL PROPERTY FORCHECK_PREPROC_SOURCES)
  get_property(_fck_preproc_sources GLOBAL PROPERTY CPP_PREPROC_SOURCES)
  #add_custom_target(all_preproc DEPENDS ${_fck_preproc_sources})
  #set_target_properties(all_preproc PROPERTIES EXCLUDE_FROM_ALL TRUE)

  #include directories for running Forcheck
  if(_fck_includes)
    string (REGEX REPLACE ";" "," _fck_incs "${_fck_includes}")
    set(_fck_incs "-I ${_fck_incs}")
  endif()
  
  # the Forcheck target
  add_custom_target(forcheck
                     COMMAND forchk -allc -rep selalib.rep -l selalib.lst ${_fck_incs} ${_fck_preproc_sources}  ${FORCHECK_EXTERNAL_FLBS}
                     WORKING_DIRECTORY  ${FORCHECK_DIRECTORY}
                     COMMENT "Running Forcheck for selalib"
                     DEPENDS all_preproc)
  set_target_properties(forcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)
endfunction()
