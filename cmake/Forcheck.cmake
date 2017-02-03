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
   set(FORCHECK_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/forcheck)
   file(MAKE_DIRECTORY ${FORCHECK_OUTPUT_DIR})

   # Some files that are required for the analysis are here
   set(FORCHECK_INPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/forcheck)

   # Forcheck library files for the external libraries
   set(FORCHECK_EXTERNAL_FLBS $ENV{FCKDIR}/share/forcheck/MPI_3.flb ${FORCHECK_INPUT_DIR}/hdf5-1_8_9.flb)
   if(NOT EXISTS "$ENV{FCKDIR}/share/forcheck/MPI_3.flb")
     message(WARNING "Forcheck: Can't find MPI_3.flb.\n Most probably the Forcheck module is not loaded.\n Try to load it and rerun cmake.")
   endif()

   MARK_AS_ADVANCED(
  FORCHECK_OUTPUT_DIR
  FORCHECK_INPUT_DIR
  FORCHECK_EXTERNAL_FLBS
  )
  
ENDIF(FORCHECK_FOUND)


# adds a custom target for running the Forcheck analysis
# call this function at the end of the CMakeLists.txt
function(add_forcheck_target)
  if(NOT FORCHECK_FOUND)
    return()
  endif()
  # retriev the lists that were creatde by add_library:
  get_property(_fck_sources GLOBAL PROPERTY CPP_SOURCES)
  get_property(_fck_includes GLOBAL PROPERTY CPP_INCLUDES)
  
  list(REMOVE_DUPLICATES _fck_sources)

  # set up include flags for preprocessing
  set(_incflags)
  foreach(i ${_fck_includes})
    set(_incflags ${_incflags} -I${i})
  endforeach()
  
  # Create custom commands for preprocessing the Fortran files
  foreach (_src ${_fck_sources})
     # Here we generate the name of the preprocessed source file
     get_filename_component(_e "${_src}" EXT)
     get_filename_component(_d "${_src}" DIRECTORY)
     get_filename_component(_n "${_src}" NAME_WE)
     string(REGEX REPLACE "F" "f" _e "${_e}")
     # get the path relative to the source dir
     string(REGEX REPLACE "^${CMAKE_SOURCE_DIR}" "" _d ${_d})
     set(_preproc_src "${CMAKE_BINARY_DIR}${_d}/${_n}_forchk${_e}")
    
     # get the compiler definitions for the file
     get_source_file_property(_defs "${_src}" COMPILE_DEFINITIONS)
     set(_defflags)
     foreach(_d ${_defs})
       set(_defflags ${_defflags} -D${_d})
     endforeach()
     
     add_custom_command(OUTPUT "${_preproc_src}"
         COMMAND gfortran  ${_incflags} ${_defflags} -cpp -E -P ${_src} | sed -f ${FORCHECK_INPUT_DIR}/sedscript -f ${FORCHECK_INPUT_DIR}/sedscript2  > ${_preproc_src}
         DEPENDS "${_src}"
         COMMENT "Preprocessing ${_src}"
         VERBATIM
       ) 
     # The preprocessed file is piped through a sed script, 
     # to break up the long lines that contain ';'.
     # To avoid trouble, first we remove comment lines that contain ';'.
     set_source_files_properties(${_preproc_src} PROPERTIES GENERATED TRUE)
     list(APPEND _fck_preproc_sources ${_preproc_src})
  endforeach()

  # group all preprocessing commands into one target
  get_property(_fck_preproc_sources GLOBAL PROPERTY CPP_PREPROC_SOURCES) #hack
  #add_custom_target(forcheck_preproc DEPENDS ${_fck_preproc_sources})
  #set_target_properties(forcheck_preproc PROPERTIES EXCLUDE_FROM_ALL TRUE)
  
  # include directories for running Forcheck
  if(_fck_includes)
    string (REGEX REPLACE ";" "," _fck_incs "${_fck_includes}")
    set(_fck_incs "-I ${_fck_incs}")
  endif()
  
  # the Forcheck target
  add_custom_target(forcheck
      COMMAND forchk -batch -allc -rep selalib.rep -l selalib.lst ${_fck_incs} ${_fck_preproc_sources}  ${FORCHECK_EXTERNAL_FLBS}
      WORKING_DIRECTORY  ${FORCHECK_OUTPUT_DIR}
      COMMENT "Running Forcheck static source code analysis"
      DEPENDS all_preproc)
  set_target_properties(forcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)

  # For each library, we will generate a forcheck target
  get_property( _library_targets GLOBAL PROPERTY LIBRARY_TARGETS )
  set(_forcheck_targets)
  foreach(_name ${_library_targets})
    get_target_property(_location ${_name} LOCATION)
    get_filename_component(_directory ${_location} DIRECTORY)
    #message(STATUS "${_directory}")
    #set(_directory $<TARGET_FILE_DIR:${_name}>)
    
    get_target_property(_sources ${_name} SOURCES)
    if (_sources)
      # we create a list of preprocessed source file names 
      set(_current_library_sources)
      foreach (_src ${_sources})
        # Here we generate the name of the preprocessed source file
        get_filename_component(_e "${_src}" EXT)
        get_filename_component(_n "${_src}" NAME_WE)
        string(REGEX REPLACE "F" "f" _e "${_e}")
        set(_preproc_src "${_directory}/${_n}${_e}")
        list(APPEND _current_library_sources ${_preproc_src})
      endforeach()
    
      # Create a list of library dependencies
      get_target_property(_link_lib ${_name} LINK_LIBRARIES)
      set(_flb_dependencies)
      foreach(_lib ${_link_lib})
        list(FIND _library_targets ${_lib} _idx)
        if(${_idx} GREATER -1)
          list(APPEND _flb_dependencies ${FORCHECK_OUTPUT_DIR}/${_lib}.flb)
        endif()
      endforeach()
      get_target_property(_link_lib2 ${_name} INTERFACE_LINK_LIBRARIES)
      if(_link_lib2)
      foreach(_lib ${_link_lib2})
        list(FIND _link_lib ${_lib} _idx)
        if(${_idx} EQUAL -1)
            message(STATUS "new LINK_LIBRARY found ${_lib}")
            list(APPEND _flb_dependencies ${FORCHECK_OUTPUT_DIR}/${_lib}.flb)
        endif()
      endforeach()
      endif()
      # add forcheck command for the library
      add_custom_command(OUTPUT ${FORCHECK_OUTPUT_DIR}/${_name}.flb
        COMMAND forchk -batch -rep ${_name}.rep -l ${_name}.lst ${_fck_incs} ${_current_library_sources} -update ${_name}.flb ${FORCHECK_EXTERNAL_FLBS} || true
        DEPENDS ${_current_library_sources} ${_flb_dependencies}
        BYPRODUCTS ${FORCHECK_OUTPUT_DIR}/${_name}.rep ${FORCHECK_OUTPUT_DIR}/${_name}.lst
        WORKING_DIRECTORY ${FORCHECK_OUTPUT_DIR}
        COMMENT Runs forcheck analysis for ${_name}
        VERBATIM
        )
    # "|| true" is needed because of Forcheck's exit status:
    # 0 no informative, warning, overflow or error messages presented
    # 2 informative, but no warning, overflow or error messages presented
    # 4 warning, but no overflow or error messages presented
    # 6 table overflow, but no error messages presented
    # 8 error messages presented
    # 16 fatal error occurred
    
      add_custom_target(forcheck_${_name} 
        DEPENDS  ${FORCHECK_OUTPUT_DIR}/${_name}.flb
      )
      set_target_properties(forcheck_${_name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
      list(APPEND _forcheck_targets ${FORCHECK_OUTPUT_DIR}/${_name}.flb)
    endif()
  endforeach()
  
  add_custom_target(forcheck_separete DEPENDS ${_forcheck_targets})
  set_target_properties(forcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)
endfunction()
