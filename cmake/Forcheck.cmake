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
  if(FORCHECK_FOUND)
    add_forcheck_target_allfiles() # process all sources at once
    add_forcheck_target_separate() # separet forcheck target for each target
  endif()
endfunction()

# Adds a single target "forcheck" to process all the source files at once.
function(add_forcheck_target_allfiles)
  # include flags needed for preprocessing the input files:
  get_property(_includes GLOBAL PROPERTY CPP_INCLUDES)
  # set up include flags for preprocessing
  set(_incflags)
  foreach(i ${_includes})
    set(_incflags ${_incflags} -I${i})
  endforeach()

  get_property(_fck_sources GLOBAL PROPERTY CPP_SOURCES)  
  list(REMOVE_DUPLICATES _fck_sources)
  
  # Create custom commands for preprocessing the Fortran files
  foreach (_src ${_fck_sources})
     # Here we generate the name of the output (preprocessed) file
     get_filename_component(_e "${_src}" EXT)
     get_filename_component(_d "${_src}" DIRECTORY)
     get_filename_component(_n "${_src}" NAME_WE)
     string(REGEX REPLACE "F" "f" _e "${_e}")
     # get the path relative to the source dir
     string(REGEX REPLACE "^${CMAKE_SOURCE_DIR}" "" _d ${_d})
     set(_preproc_name "${CMAKE_BINARY_DIR}${_d}/${_n}_forchk${_e}")
    
     # get the compiler definitions for the file
     get_source_file_property(_defs "${_src}" COMPILE_DEFINITIONS)
     set(_defflags)
     foreach(_d ${_defs})
       set(_defflags ${_defflags} -D${_d})
     endforeach()
     
     # Create the preprocessor command
     # The preprocessed file is piped through a sed script, 
     # to break up the long lines that contain ';'.
     add_custom_command(OUTPUT "${_preproc_name}"
         COMMAND gfortran  ${_incflags} ${_defflags} -cpp -E -P ${_src} | sed -f ${FORCHECK_INPUT_DIR}/sedscript -f ${FORCHECK_INPUT_DIR}/sedscript2  > ${_preproc_name}
         DEPENDS "${_src}"
         COMMENT "Preprocessing ${_src}"
         VERBATIM
       ) 
     list(APPEND _fck_preproc_sources ${_preproc_name})
  endforeach()

  # group all preprocessing commands into one target
  get_property(_fck_preproc_sources GLOBAL PROPERTY CPP_PREPROC_SOURCES) #hack
  #add_custom_target(forcheck_preproc DEPENDS ${_fck_preproc_sources})
  #set_target_properties(forcheck_preproc PROPERTIES EXCLUDE_FROM_ALL TRUE)
  
  get_forcheck_includes(_fck_incflags)
  # the Forcheck target
  add_custom_target(forcheck
      COMMAND forchk -batch -allc -rep selalib.rep -l selalib.lst ${_fck_incflags} ${_fck_preproc_sources}  ${FORCHECK_EXTERNAL_FLBS}
      WORKING_DIRECTORY  ${FORCHECK_OUTPUT_DIR}
      COMMENT "Running Forcheck static source code analysis"
      DEPENDS all_preproc)
  set_target_properties(forcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)
endfunction()

# Create separate Forcheck target for each compilation target
# The forcheck targets are named forcheck_NAME, where NAME is
# the name of any target added by ADD_LIBRARY or ADD_EXECUTABLE
# Additionally, single target is added to check all libraries
# individually: forcheck_separate
function(add_forcheck_target_separate)
  get_forcheck_includes(_fck_incflags)
  get_property(_library_targets GLOBAL PROPERTY LIBRARY_TARGETS)
  set(_forcheck_targets)
  foreach(_name ${_library_targets})
    get_target_property(_directory ${_name} SOURCE_DIR)
    #get_target_property(_location ${_name} LOCATION)
    #get_filename_component(_directory ${_location} DIRECTORY)
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
      get_flb_dependencies(${_name} _flb_dependencies)
      string(REPLACE ";" " " _flbs "${_flb_dependencies}")
      # add forcheck command for the library
      add_custom_command(OUTPUT ${FORCHECK_OUTPUT_DIR}/${_name}.flb
        COMMAND forchk -allc -batch -rep ${_name}.rep -l ${_name}.lst ${_fck_incflags} ${_current_library_sources} -update ${_name}.flb ${_flbs} ${FORCHECK_EXTERNAL_FLBS} || true
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
    
      add_custom_target(forcheck_${_name} DEPENDS  ${FORCHECK_OUTPUT_DIR}/${_name}.flb)
      set_target_properties(forcheck_${_name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
      list(APPEND _forcheck_targets ${FORCHECK_OUTPUT_DIR}/${_name}.flb)
    endif()
  endforeach()
  
  add_custom_target(forcheck_separate DEPENDS ${_forcheck_targets})
  set_target_properties(forcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)
endfunction()

# Returns a comma separated list of all the preprocessor flags
function(get_forcheck_includes _output_name)
  get_property(_fck_includes GLOBAL PROPERTY CPP_INCLUDES)
  if(_fck_includes)
    string (REGEX REPLACE ";" "," _fck_incs "${_fck_includes}")
    set(_fck_incs "-I ${_fck_incs}")
    set(${_output_name} ${_fck_incs} PARENT_SCOPE)
  endif()
endfunction()

# Get a list of all library dependencies (including transitive dependecies).
# It won't work properly if we have generator expressions in LINK_LIBRARIES.
function(get_flb_dependencies _name _output_name)   
  # we will cross check the library names with the LIBRARY_TARGETS list
  get_property( _library_targets GLOBAL PROPERTY LIBRARY_TARGETS)
    
  set(_dependencies) # the output will be stored here
  # we do a breadth first search, and list all the LINK_LIBRARIES that
  # are stored in the following two properties:
  get_target_property(_link_lib ${_name} LINK_LIBRARIES)
  get_target_property(_iflink_lib ${_name} INTERFACE_LINK_LIBRARIES)    
  
  list(APPEND _link_lib ${_iflink_lib})
  list(LENGTH _link_lib _len)
  set(_idx 0)
  while(_idx LESS _len)
    list(GET _link_lib ${_idx} _libname)
    # check if it is a library target that is part of selalib
    list(FIND _library_targets ${_libname} _tmpidx)
    if(${_tmpidx} GREATER -1)
      list(FIND _dependencies ${_libname} _tmpidx)
      if(${_tmpidx} EQUAL -1)
        # not yet in the output list, we include it
        list(APPEND _dependencies ${_libname})
        # now extend the _link_lib list with dependencies of _lib
        get_target_property(_tmp ${_libname} LINK_LIBRARIES)
        get_target_property(_tmp2 ${_libname} INTERFACE_LINK_LIBRARIES)
        list(APPEND _tmp ${_tmp2})
        foreach(_libname2 ${_tmp})
          # we include only those that are not yet in _link_lib
          list(FIND _link_lib ${_libname2} _tmpidx)
          if(${_tmpidx} EQUAL -1)
            list(APPEND _link_lib ${_libname2})
          endif()
        endforeach()
      endif()
    endif()
    list(LENGTH _link_lib _len)
    math(EXPR _idx "${_idx} + 1")
  endwhile()
  
  #message(STATUS "${_name} depends on ${_dependencies}")
  set(_flb_dependencies)
  foreach(_libname ${_dependencies})
    list(APPEND _flb_dependencies ${FORCHECK_OUTPUT_DIR}/${_libname}.flb)
  endforeach()
  set(${_output_name} ${_flb_dependencies} PARENT_SCOPE)
endfunction()