# CMake module to create targets to run Forcheck static source code analysis
#
# Usage: Add the following line to the beginning of the main CMakeLists.txt file
# include(PreprocessorTarget)
#
# At the end of the same CMakeLists.txt call add_forcheck_target()
#
# After generating the makefiles one can use make forcheck
#
# This will run a Forcheck analysis over all the source files, and the results
# of the analysis will be saved forcheck/selalib.lst. Additionally
# forcheck/selalib.rep file will be also created (it does not contain anything
# new compared to the lst file, the rep file is just a shorter report without
# source listing and without message summary).
#
# If you call cmake with -DFORCHECK_SEPARATE_TARGETS=1, then for each target, an
# additional forcheck target will be also created with the name
# forcheck_TARGETNAME. For example to run Forcheck for sll_working_precision use
#
# make forcheck_sll_working_precision
#
# The results will be saved in forcheck/sll_working_precision.lst Additionally
# forcheck/sll_working_precision.flb will be also generated which stores
# interface information about the library. It is used to facilitate testing
# other libraries which depend on sll_working_precision.
#
# To run forcheck for all library individually with a single command use make
# forcheck_separate

find_program(
  FORCHECK_EXECUTABLE
  NAMES forchk
  PATHS $ENV{FCKDIR}/bin /usr/bin /bin /usr/local/bin
  DOC "Performs a full static analysis of Fortran programs.")

if(FORCHECK_EXECUTABLE)
  set(FORCHECK_FOUND "YES")
endif(FORCHECK_EXECUTABLE)

mark_as_advanced(FORCHECK_FOUND FORCHECK_EXECUTABLE)

if(FORCHECK_FOUND)
  # The result of the forcheck analysis will go into this directory
  set(FORCHECK_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/forcheck)
  file(MAKE_DIRECTORY ${FORCHECK_OUTPUT_DIR})

  # Some files that are required for the analysis are here
  set(FORCHECK_INPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/forcheck)

  # Forcheck library files for the external libraries
  set(FORCHECK_EXTERNAL_FLBS $ENV{FCKDIR}/share/forcheck/MPI_3.flb
                             ${FORCHECK_INPUT_DIR}/hdf5-1_8_18.flb)
  if(NOT EXISTS "$ENV{FCKDIR}/share/forcheck/MPI_3.flb")
    message(
      WARNING
        "Forcheck: Can't find MPI_3.flb.\n Most probably the Forcheck module is not loaded.\n Try to load it and rerun cmake."
    )
  endif()

  # Trying to set Fortran compiler emulation, somehow it does not work
  set(ENV{FCKCNF} "$ENV{FCKDIR}/share/forcheck/intel15.cnf") # Fortran 2003 with
                                                             # Intel extensions
  # Alternatively: set (ENV{FCKCNF} "$ENV{FCKDIR}/share/forcheck/f03.cnf") #
  # Fortran 2003 without extensions set (ENV{FCKCNF}
  # "$ENV{FCKDIR}/share/forcheck/f08.cnf") # Fortran 2008
  mark_as_advanced(FORCHECK_OUTPUT_DIR FORCHECK_INPUT_DIR
                   FORCHECK_EXTERNAL_FLBS)

endif(FORCHECK_FOUND)

# Adds custom targets for running the Forcheck analysis call this function at
# the end of the CMakeLists.txt
function(add_forcheck_target)
  if(FORCHECK_FOUND)
    add_forcheck_target_allfiles() # process all sources at once
    if(FORCHECK_SEPARATE_TARGETS)
      add_forcheck_target_separate()
      # separet forcheck target for each library/executable
    endif()
  endif()
endfunction()

# Adds a single target "forcheck" to process all the source files at once.
function(add_forcheck_target_allfiles)
  # Forcheck will analyze already preprocessed files. These files will not
  # contain #include proprocessor directives anymore, but they can still have
  # fortran include lines. Therefore, we need to get the include flags and pass
  # it to Forcheck
  get_forcheck_includes(_fck_incflags)

  # List of preprocessed source files
  get_property(_fck_preproc_sources GLOBAL PROPERTY CPP_PREPROC_SOURCES)

  # The Forcheck target
  add_custom_target(
    forcheck
    COMMAND forchk -batch -allc -rep selalib.rep -l selalib.lst ${_fck_incflags}
            ${_fck_preproc_sources} ${FORCHECK_EXTERNAL_FLBS} || true
    WORKING_DIRECTORY ${FORCHECK_OUTPUT_DIR}
    COMMENT "Running Forcheck static source code analysis"
    DEPENDS ${_fck_preproc_sources})
  # "|| true" is used because of Forcheck's exit status: 0 no informative,
  # warning, overflow or error messages presented 2 informative, but no warning,
  # overflow or error messages presented 4 warning, but no overflow or error
  # messages presented 6 table overflow, but no error messages presented 8 error
  # messages presented 16 fatal error occurred
  set_target_properties(forcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)
endfunction()

# Create separate Forcheck target for each compilation target The forcheck
# targets are named forcheck_NAME, where NAME is the name of any target added by
# ADD_LIBRARY or ADD_EXECUTABLE Additionally, single target is added to check
# all libraries individually: forcheck_separate
function(add_forcheck_target_separate)
  get_forcheck_includes(_fck_incflags)

  get_property(_target_list GLOBAL PROPERTY LIBRARY_TARGETS)
  # get_property(_executable_list GLOBAL PROPERTY EXECUTABLE_TARGETS)
  # list(APPEND _target_list ${_executable_list})
  set(_forcheck_targets)
  # todo: append the list of executables to target_list
  foreach(_name ${_target_list})
    get_target_property(_source_directory ${_name} SOURCE_DIR)
    get_target_property(_sources ${_name} SOURCES)
    if(_sources)
      list(REMOVE_DUPLICATES _sources)
      # we create a list of preprocessed source file names
      set(_current_library_sources)
      foreach(_src ${_sources})
        get_preprocessed_filename(${_source_directory}/${_src} _preproc_src)
        list(APPEND _current_library_sources ${_preproc_src})
      endforeach()

      # Create a list of library dependencies
      get_flb_dependencies(${_name} _flb_dependencies)
      # The COMMAND option will need a space separated list while the DEPNDS
      # option nedds a comma separated list
      string(REPLACE ";" " " _flb_space "${_flb_dependencies}")

      # Add forcheck command for the library
      add_custom_command(
        OUTPUT ${FORCHECK_OUTPUT_DIR}/${_name}.flb
        COMMAND
          forchk -allc -batch -rep ${_name}.rep -l ${_name}.lst ${_fck_incflags}
          ${_current_library_sources} -update ${_name}.flb ${_flb_space}
          ${FORCHECK_EXTERNAL_FLBS} || true
        DEPENDS ${_current_library_sources} ${_flb_dependencies}
        BYPRODUCTS ${FORCHECK_OUTPUT_DIR}/${_name}.rep
                   ${FORCHECK_OUTPUT_DIR}/${_name}.lst
        WORKING_DIRECTORY ${FORCHECK_OUTPUT_DIR}
        COMMENT Runs forcheck analysis for ${_name}
        VERBATIM)
      # "|| true" is needed because of Forcheck's exit status: 0 no informative,
      # warning, overflow or error messages presented 2 informative, but no
      # warning, overflow or error messages presented 4 warning, but no overflow
      # or error messages presented 6 table overflow, but no error messages
      # presented 8 error messages presented 16 fatal error occurred
      add_custom_target(forcheck_${_name}
                        DEPENDS ${FORCHECK_OUTPUT_DIR}/${_name}.flb)
      set_target_properties(forcheck_${_name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
      list(APPEND _forcheck_targets ${FORCHECK_OUTPUT_DIR}/${_name}.flb)
    endif()
  endforeach()

  add_custom_target(forcheck_separate DEPENDS ${_forcheck_targets})
  set_target_properties(forcheck PROPERTIES EXCLUDE_FROM_ALL TRUE)
endfunction()

# Returns a list of all the include flags used by any target. The list is in a
# format that can be given to Forcheck
function(get_forcheck_includes _output_name)
  # Get the list of all targes
  get_property(_target_list GLOBAL PROPERTY LIBRARY_TARGETS)
  get_property(_executable_list GLOBAL PROPERTY EXECUTABLE_TARGETS)
  list(APPEND _target_list ${_executable_list})

  # Collect all the include directories
  set(_includes)
  foreach(_name ${_target_list})
    get_target_property(_dirs ${_name} INCLUDE_DIRECTORIES)
    if(_dirs)
      list(APPEND _includes "${_dirs}")
    endif()
    get_target_property(_dirs ${_name} INTERFACE_INCLUDE_DIRECTORIES)
    if(_dirs)
      list(APPEND _includes "${_dirs}")
    endif()
    list(REMOVE_DUPLICATES _includes)
  endforeach()

  # Transform it to a Format that forcheck accepts
  if(_includes)
    string(REGEX REPLACE ";" "," _fck_incs "${_includes}")
    set(_fck_incs "-I ${_fck_incs}")
    set(${_output_name}
        ${_fck_incs}
        PARENT_SCOPE)
  endif()
endfunction()

# Get a list of all library dependencies (including transitive dependecies). It
# won't work properly if we have generator expressions in LINK_LIBRARIES.
function(get_flb_dependencies _name _output_name)
  # we will cross check the library names with the LIBRARY_TARGETS list
  get_property(_library_targets GLOBAL PROPERTY LIBRARY_TARGETS)

  set(_dependencies) # the output will be stored here
  # we do a breadth first search, and list all the LINK_LIBRARIES that are
  # stored in the following two properties:
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

  set(_flb_dependencies)
  foreach(_libname ${_dependencies})
    list(APPEND _flb_dependencies ${FORCHECK_OUTPUT_DIR}/${_libname}.flb)
  endforeach()
  set(${_output_name}
      ${_flb_dependencies}
      PARENT_SCOPE)
endfunction()
