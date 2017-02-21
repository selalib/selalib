# CMake module to create a target to preprocess all the fortran libraries
#
# Usage:
# Add the following line to the beginning of the main CMakeLists.txt file
# include(PreprocessorTarget)
# 
# At the end of the same CMakeLists.txt call
# add_preprocessor_target()
# 
# After generating the makefiles one can call
# make all_preproc
#
# This module overrides the add_library cmake command to create a list of source
# files that can be given to Forcheck for analysis. Therefore, we should place the
# include(PreprocessorTarget) command before any add_library or add_subdirectory commands.
#
# The source files are preprocessed individually using the C preprocessor.
# The add_preprocessor_target() function generates the commands for the
# preprocessor. It should be called after all the libraries and 
# subdirectories are included.
# 
# Author of ForcheckTargets.cmake: Tamas Feher <tamas.bela.feher@ipp.mpg.de>
#
# Modifications
# -------------
#   - 22 Oct 2015: only preprocess files (Yaman Güçlü [YG], IPP Garching).
#   - 02 Nov 2015: also preprocess executable sources (YG).
#   - 26 Nov 2015: add OpenMP flag (YG).
#   - 02 Dec 2015: fix dependency bug (YG).
#   - 15 Jan 2016: store names of all libraries (YG).
#   - 19 Jan 2016: 'collect_source_info' handles libraries with no sources (YG)
#   - 03 Jan 2017: Add PGI compiler (PN)

if(__add_all_preproc)
   return()
endif()

set(__add_all_preproc YES)

# List of targets created by "add_library" instructions
set_property(GLOBAL PROPERTY LIBRARY_TARGETS "")
set_property(GLOBAL PROPERTY EXECUTABLE_TARGETS "")

# List of source files to be analyzed
#set_property(GLOBAL PROPERTY CPP_SOURCES "")
set_property(GLOBAL PROPERTY CPP_PREPROC_SOURCES "")

# List of include directories
#set_property(GLOBAL PROPERTY CPP_INCLUDES "")

# Preprocessor flags
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(preprocessor_only_flags -EP)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(preprocessor_only_flags -cpp -E -P)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  set(preprocessor_only_flags -Mcpp -E -P )
else()
  message(SEND_ERROR "Unknown preprocessor flags for current compiler")
endif()

#==============================================================================
# FUNCTION: collect_source_info
#==============================================================================
# Create a list of source files, to be later used to run the preprocessor
function( collect_source_info _name )

endfunction()

#==============================================================================
# FUNCTION: collect_library_name
#==============================================================================
# Collect names of all library targets.
function( collect_library_name _name )
  get_property( _library_targets GLOBAL PROPERTY LIBRARY_TARGETS )
  list( APPEND _library_targets "${_name}" )
  set_property( GLOBAL PROPERTY LIBRARY_TARGETS ${_library_targets} )
endfunction()

#==============================================================================
# FUNCTION: collect_executable_name
#==============================================================================
# Collect names of all library targets.
function( collect_executable_name _name )
  get_property( _targets GLOBAL PROPERTY EXECUTABLE_TARGETS )
  list( APPEND _targets "${_name}" )
  set_property( GLOBAL PROPERTY EXECUTABLE_TARGETS ${_targets} )
endfunction()

#==============================================================================
# FUNCTION: store_current_dir
#==============================================================================
# Add property to target: directory with the currently processed CMakeLists.txt
function( store_current_dir _name )
  set_target_properties( ${_name} PROPERTIES
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
endfunction()

#==============================================================================
# FUNCTION: add_library
#==============================================================================
# We override the add_library built in function.
function( add_library _name )
  _add_library( ${_name} ${ARGN} ) # Call the original function
  collect_library_name( ${_name} ) # Store library name in proper list
  collect_source_info ( ${_name} ) # Create a list of source files
  store_current_dir   ( ${_name} ) # Store current directory in target property
endfunction()

#==============================================================================
# FUNCTION: add_executable
#==============================================================================
# We override the add_executable built in function.
function( add_executable _name )
  _add_executable( ${_name} ${ARGN} ) # Call the original function
  collect_executable_name( ${_name} )
  collect_source_info( ${_name} )     # Create a list of source files
  store_current_dir ( ${_name} )
endfunction()

#==============================================================================
# FUNCTION: add_preprocessor_target
#==============================================================================
# adds a custom target for running the C preprocessor on all source files
# call this function at the end of the CMakeLists.txt
function(add_preprocessor_target)

  # If needed, add OpenMP flag to preprocessor flags
  if(OPENMP_ENABLED)
    set(preprocessor_only_flags ${preprocessor_only_flags} ${OpenMP_Fortran_FLAGS})
  endif()

  get_property(_target_list GLOBAL PROPERTY LIBRARY_TARGETS)
  get_property(_executable_list GLOBAL PROPERTY EXECUTABLE_TARGETS)
  list(APPEND _target_list ${executable_list})
  
  set(_preproc_sources) # we collect the names of the preprocessed source files here
  
  foreach(_name ${_target_list})
    message(STATUS "checking source files for target ${_name}")
    get_target_property(_source_directory ${_name} SOURCE_DIR)
    get_target_property(_target_defs ${_name} COMPILE_DEFINITIONS)
    get_property(_dir_defs  DIRECTORY ${_source_directory} PROPERTY COMPILE_DEFINITIONS)
    if(_target_defs)
      message(STATUS "target defs for ${_name} ${_dir_defs}")
    else()
    set(_target_defs)
    endif()
    if(_dir_defs)
      list(APPEND _target_defs ${_dir_defs})
      message(STATUS "dir defs for ${_name} ${_dir_defs}")
    endif()
    get_include_flags(${_name} _incflags) # set up include flags specific for each target
    
    get_target_property(_sources ${_name} SOURCES)
    # Problems with the preprocessed files:
    #  - they can contain arbitrary long lines because of macro expansion,
    #  - Intel's proprocessor brakes up these lines, but create non-standard
    #    conforming code.
    # Therefore in the following loop we create custom commands to preprocess the 
    # source files with gfortran and fix the problem with long lines using sed.
  
    # now we create targets for preprocessing the files
    if(_sources)
     foreach (_src ${_sources})
      check_if_fortran_file(${_src} _fortran_file)
      # get_source_file_property(_lang "${_src_loc}" LANGUAGE) does not work because we are not in 
      # same directory where the target was added, therefore we can only see the LOCATION property
      if(_fortran_file)   
        set(_src_loc "${_source_directory}/${_src}")
        get_preprocessed_filename(${_src_loc} _preproc_name)
        list(FIND _preproc_sources ${_preproc_name} _tmpidx)
        if (${_tmpidx} EQUAL -1)
          # not yet in the list of preprocessed files create a target
          get_compile_definitions(${_src_loc} "${_target_defs}" _defflags)
          # Create the preprocessor command
          # The preprocessed file is piped through a sed script, 
          # to break up the long lines that contain ';'.
          # To avoid trouble, we delete comment lines that contain  ';'.
          add_custom_command(OUTPUT "${_preproc_name}"
            COMMAND gfortran  ${_incflags} ${_defflags} -cpp -E -P ${_src_loc} | sed -e "/^.\\{132\\}/s/!.*/ /" -e "/^.\\{132\\}/s/; */\\n/g" > ${_preproc_name}
            #COMMAND ${CMAKE_Fortran_COMPILER} ${incflags} ${defflags} ${preprocessor_only_flags} ${_src} > ${_preproc_src}
            DEPENDS "${_src_loc}"
            COMMENT "Preprocessing ${_src}"
            VERBATIM
          ) 
          list(APPEND _preproc_sources ${_preproc_name})
        else()
          message(STATUS "file ${_src_loc} is already processed with idx ${_tmpidx}")
        endif()
      else()
        message(STATUS "${_src} is not a fortran file")
      endif()
     endforeach()
    endif()
  endforeach()

  # Group all preprocessing commands into one target
  set_property(GLOBAL PROPERTY CPP_PREPROC_SOURCES ${_preproc_sources})
  add_custom_target(all_preproc DEPENDS ${_preproc_sources})
  set_target_properties(all_preproc PROPERTIES EXCLUDE_FROM_ALL TRUE)

  # Clean up *.s files just generated (if any)
  add_custom_command( TARGET all_preproc POST_BUILD
    COMMAND find . -iname "*.s" -type f -delete
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ) 
endfunction()

# checks the extension of the file _src to determine whether it is a fortran file
function(check_if_fortran_file _src _output_name)
  get_filename_component(_ext "${_src}" EXT)
  string(TOLOWER ${_ext} _ext)
  set(_fortran_file FALSE)
  if ("${_ext}" STREQUAL ".f")
    set(_fortran_file TRUE)
  elseif("${_ext}" STREQUAL ".f90")
    set(_fortran_file TRUE)
  endif()
  set(${_output_name} ${_fortran_file} PARENT_SCOPE)
endfunction()
# get the list of include directories for _target, and transorm into
# a space separated list with -I prefix, to be given as proprocessor flags
function(get_include_flags _target _output_name)
    set(_includes)
    get_target_property(_dirs ${_target} INCLUDE_DIRECTORIES)
    if(_dirs)
      list(APPEND _includes "${_dirs}")
    endif()
    get_target_property(_dirs ${_target} INTERFACE_INCLUDE_DIRECTORIES)
    if(_dirs)
      list(APPEND _includes "${_dirs}")
    endif()
    if(_includes)
      list(REMOVE_DUPLICATES _includes)
    endif()
    set(_incflags)
    foreach(i ${_includes})
      set(_incflags ${_incflags} -I${i})
    endforeach()
  set(${_output_name} "${_incflags}" PARENT_SCOPE)
endfunction()

# get a list of compile definition prefixed with -D, separated by space
function(get_compile_definitions _src_loc _target_defs _output_name)
  set(_defs)
  if(_target_defs)
    list(APPEND _defs ${_target_defs})
  endif()
  get_source_file_property(_cpp_defs "${_src_loc}" COMPILE_DEFINITIONS)
  #unfortunatelly COMPILE_DEFINITIONS are not visible outside the directory 
  #where the target was added
  if(_cpp_defs)
    list(APPEND _defs "${_cpp_defs}")
  endif()
  if(_defs)
    list(REMOVE_DUPLICATES _defs)
  endif()
  message(STATUS "${_src_loc} ${_defs}")
  set(_defflags)
  foreach(_d ${_defs})
    set(_defflags ${_defflags} -D${_d})
  endforeach()
  set(${_output_name} "${_defflags}" PARENT_SCOPE)
endfunction() 