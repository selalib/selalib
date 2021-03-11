# CMake module to create a target to preprocess all the fortran libraries
#
# Usage: Add the following line to the beginning of the main CMakeLists.txt file
# include(PreprocessorTarget)
#
# At the end of the same CMakeLists.txt call add_preprocessor_target()
#
# After generating the makefiles one can call make all_preproc
#
# This module overrides the add_library and add_executable cmake commands to
# create a list of targets. The list of targets is later used to preprocess
# their source files. We should place the include(PreprocessorTarget) command
# before any add_library or add_subdirectory commands.
#
# The source files are preprocessed individually using the Fortran compiler.
#
# The add_preprocessor_target() function generates the commands for the
# preprocessor. It should be called after all the libraries and subdirectories
# are included.
#
# Author of ForcheckTargets.cmake: Tamas Feher <tamas.bela.feher@ipp.mpg.de>
#
# Modifications
# -------------
# * 22 Oct 2015: only preprocess files (Yaman Güçlü [YG], IPP Garching).
# * 02 Nov 2015: also preprocess executable sources (YG).
# * 26 Nov 2015: add OpenMP flag (YG).
# * 02 Dec 2015: fix dependency bug (YG).
# * 15 Jan 2016: store names of all libraries (YG).
# * 19 Jan 2016: 'collect_source_info' handles libraries with no sources (YG)
# * 03 Jan 2017: Add PGI compiler (PN)
# * 22 Feb 2017: To be able to use the preprocessed files for Forcheck analysis,
#   the following changes were implemented: - Store a list of executables too. -
#   No need to store list of source files, it can be retreived from the lists of
#   targets (by combining the SOURCES and the SOURCE_DIR property). -
#   Collect_source_info not needed anymore. - We set the include flags and
#   compile  definitions for each target individually. - Preprocessed files can
#   be piped through a sed script to break long lines. - If Forcheck is
#   available, do not use Intel compiler to preprocess the source files. The
#   Intel preprocessor brakes long lines in non standard conforming way. We try
#   GFortran instead. - For older versions of GFortran, only copy the .f[90]
#   files. (Tamas Feher)

if(__add_all_preproc)
  return()
endif()

set(__add_all_preproc YES)

# List of targets created by "add_library"  and "add_executable" instructions
set_property(GLOBAL PROPERTY LIBRARY_TARGETS "")
set_property(GLOBAL PROPERTY EXECUTABLE_TARGETS "")

# List of source files created during preprocessing
set_property(GLOBAL PROPERTY CPP_PREPROC_SOURCES "")

# Whether to proprocess or copy .f or .f90 files (.F and .F90 will be always
# preprocessed)
set(copy_lowercase_f_files FALSE)

# Whether we should brake too long lines in the preprocessed files (needed for
# Forcheck)
set(break_long_lines FALSE)

# Preprocessor command and flags
set(preprocessor_command ${CMAKE_Fortran_COMPILER})
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(preprocessor_only_flags -EP)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(preprocessor_only_flags -cpp -E -P)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
  set(preprocessor_only_flags -Mcpp -E -P)
else()
  message(SEND_ERROR "Unknown preprocessor flags for current compiler")
endif()

if(FORCHECK_FOUND)
  set(break_long_lines TRUE)
  if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    if(${Fortran_COMPILER_VERSION} VERSION_LESS 4.9)
      # Earlier versions of GFortran did not process the .f files, so we will
      # copy these instead of preprocessing.
      set(copy_lowercase_f_files TRUE)
    endif()
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES Intel)
    # Override proprocessor command because Intel creates non standard
    # conforming code (it breaks long lines incorrectly).
    find_program(GFORTRAN_COMPILER gfortran
                 DOC "GFortran is needed for preprocessing")
    if(GFORTRAN_COMPILER)
      set(preprocessor_command ${GFORTRAN_COMPILER})
      set(preprocessor_only_flags -cpp -E -P)
      message(STATUS "Preprocessor command changed to ${GFORTRAN_COMPILER}")
      # check the version
      execute_process(COMMAND ${GFORTRAN_COMPILER} --version
                      OUTPUT_VARIABLE _tmp)
      string(REGEX MATCH "[4-7]\\.[0-9]\\.[0-9]" GFortran_COMPILER_VERSION
                   ${_tmp})
      if(${GFortran_COMPILER_VERSION} VERSION_LESS 4.9)
        set(copy_lowercase_f_files TRUE)
      endif()
    else()
      message(
        WARNING
          "Forcheck analysis will have problem with ifort preprocessed files")
    endif()
  endif()
endif()

# ==============================================================================
# FUNCTION: collect_library_name
# ==============================================================================
# Collect names of all library targets.
function(collect_library_name _name)
  get_property(_library_targets GLOBAL PROPERTY LIBRARY_TARGETS)
  list(APPEND _library_targets "${_name}")
  set_property(GLOBAL PROPERTY LIBRARY_TARGETS ${_library_targets})
endfunction()

# ==============================================================================
# FUNCTION: collect_executable_name
# ==============================================================================
# Collect names of all library targets.
function(collect_executable_name _name)
  get_property(_targets GLOBAL PROPERTY EXECUTABLE_TARGETS)
  list(APPEND _targets "${_name}")
  set_property(GLOBAL PROPERTY EXECUTABLE_TARGETS ${_targets})
endfunction()

# ==============================================================================
# FUNCTION: store_current_dir
# ==============================================================================
# Add property to target: directory with the currently processed CMakeLists.txt
function(store_current_dir _name)
  set_target_properties(${_name} PROPERTIES SOURCE_DIR
                                            ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

# ==============================================================================
# FUNCTION: add_library
# ==============================================================================
# We override the add_library built in function.
function(add_library _name)
  _add_library(${_name} ${ARGN}) # Call the original function
  collect_library_name(${_name}) # Store library name in proper list
  store_current_dir(${_name}) # Store current directory in target property
endfunction()

# ==============================================================================
# FUNCTION: add_executable
# ==============================================================================
# We override the add_executable built in function.
function(add_executable _name)
  _add_executable(${_name} ${ARGN}) # Call the original function
  collect_executable_name(${_name})
  store_current_dir(${_name})
endfunction()

# ==============================================================================
# FUNCTION: add_preprocessor_target
# ==============================================================================
# adds a custom target for running the C preprocessor on all source files call
# this function at the end of the CMakeLists.txt
function(add_preprocessor_target)
  # If needed, add OpenMP flag to preprocessor flags
  if(OPENMP_ENABLED AND ${preprocessor_command} STREQUAL
                        ${CMAKE_Fortran_COMPILER})
    set(preprocessor_only_flags ${preprocessor_only_flags}
                                ${OpenMP_Fortran_FLAGS})
    # else()
    # we have modified the preprocessor command because of Forcheck. We do not
    # add the OpenMP flag, since it might not work with the modified
    # preprocessor command
  endif()

  get_property(_target_list GLOBAL PROPERTY LIBRARY_TARGETS)
  get_property(_executable_list GLOBAL PROPERTY EXECUTABLE_TARGETS)
  list(APPEND _target_list ${_executable_list})

  # set(_all_sources)   # the list of all source files, not needed now
  set(_preproc_sources) # we collect the names of the preprocessed source files

  foreach(_name ${_target_list})
    # set up include flags and compile definitions specific for each target
    get_compile_definitions(${_name} _defflags)
    get_include_flags(${_name} _incflags)
    get_target_property(_source_directory ${_name} SOURCE_DIR)
    get_target_property(_sources ${_name} SOURCES)
    # now we create targets for preprocessing the files
    if(_sources)
      foreach(_src ${_sources})
        check_source_file(${_src} _fortran_file _lowercase_f)
        if(_fortran_file)
          set(_src_loc "${_source_directory}/${_src}")
          get_preprocessed_filename(${_src_loc} _preproc_name)
          list(FIND _preproc_sources ${_preproc_name} _tmpidx)
          if(${_tmpidx} EQUAL -1)
            # file was not yet processed list(APPEND _all_sources ${_src_loc})
            list(APPEND _preproc_sources ${_preproc_name})

            if(${copy_lowercase_f_files} AND ${_lowercase_f})
              # just copy the file
              add_custom_command(
                OUTPUT ${_preproc_name}
                COMMAND ${CMAKE_COMMAND} -E copy ${_src_loc} ${_preproc_name}
                DEPENDS "${_src_loc}"
                COMMENT "Copying ${_src} (no need to preprocess)")
            elseif(break_long_lines)
              # The preprocessed file can contain arbitrary long lines because
              # of macro expansion. We fix the problem using sed.
              #
              # sed -e "/^.\\{132\\}/s/!.*/ /" matches all lines that are at
              # least 132 characters long, and replaces everything afte '!'
              # character with space. The intention is to delete comments from
              # long lines, because braking comment lines would cause trouble.
              # This sed command would not work for long lines that have a
              # character context with '!', but currently SeLaLib has no such
              # lines.
              #
              # sed -e "/^.\\{132\\}/s/; */\\n/g" matches all lines that are at
              # least 132 characters long, and replaces every ';' character with
              # newline.

              # Create the preprocessor command
              add_custom_command(
                OUTPUT "${_preproc_name}"
                COMMAND
                  ${preprocessor_command} ${_incflags} ${_defflags}
                  ${preprocessor_only_flags} ${_src_loc} | sed -e
                  "/^.\\{132\\}/s/!.*/ /" -e "/^.\\{132\\}/s/; */\\n/g" >
                  ${_preproc_name}
                DEPENDS "${_src_loc}"
                COMMENT "Preprocessing ${_src}"
                VERBATIM)
            else()
              # just preprocess without using any sed script
              add_custom_command(
                OUTPUT "${_preproc_name}"
                COMMAND
                  ${preprocessor_command} ${_incflags} ${_defflags}
                  ${preprocessor_only_flags} ${_src_loc} > ${_preproc_name}
                DEPENDS "${_src_loc}"
                COMMENT "Preprocessing ${_src}"
                VERBATIM)
            endif()
          endif()
        endif()
      endforeach()
    endif()
  endforeach()

  # Group all preprocessing commands into one target set_property(GLOBAL
  # PROPERTY CPP_SOURCES ${_all_sources})
  set_property(GLOBAL PROPERTY CPP_PREPROC_SOURCES ${_preproc_sources})
  add_custom_target(all_preproc DEPENDS ${_preproc_sources})
  set_target_properties(all_preproc PROPERTIES EXCLUDE_FROM_ALL TRUE)

  # Clean up *.s files just generated (if any)
  add_custom_command(
    TARGET all_preproc
    POST_BUILD
    COMMAND find . -iname "*.s" -type f -delete
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endfunction()

# ==============================================================================
# FUNCTION: check_if_fortran_file
# ==============================================================================
# Checks the extension of the file to determine whether it is a fortran file.
# Girst argument will be TRUE if _src is a fortran file (.f, .F, .f90, or .F90),
# second argument will be TRUE if extension is .f or .f90, otherwise returns
# FALSE.
function(check_source_file _src _output_name1 _output_name2)
  # get_source_file_property(_lang "${_src_loc}" LANGUAGE) does not work because
  # we are not in same directory where the target was added (we can only see the
  # LOCATION property).
  get_filename_component(_ext "${_src}" EXT)
  set(_fortran_file FALSE)
  set(_lowercase_f FALSE)
  if("${_ext}" STREQUAL ".f" OR "${_ext}" STREQUAL ".f90")
    set(_fortran_file TRUE)
    set(_lowercase_f TRUE)
  elseif("${_ext}" STREQUAL ".F" OR "${_ext}" STREQUAL ".F90")
    set(_fortran_file TRUE)
  endif()
  set(${_output_name1}
      ${_fortran_file}
      PARENT_SCOPE)
  set(${_output_name2}
      ${_lowercase_f}
      PARENT_SCOPE)
endfunction()

# ==============================================================================
# FUNCTION: get_include_flags
# ==============================================================================
# get the list of include directories for _target, and transorm it into a space
# separated list with -I prefix, to be given as proprocessor flags
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
  set(${_output_name}
      "${_incflags}"
      PARENT_SCOPE)
endfunction()

# ==============================================================================
# FUNCTION: get_compile_definitions
# ==============================================================================
# get a list of compile definition for a target prefixed with -D, separated by
# space
function(get_compile_definitions _name _output_name)
  get_target_property(_target_defs ${_name} COMPILE_DEFINITIONS)
  if(NOT _target_defs)
    set(_target_defs) # to avoid creating -DNOTFOUND flags
  endif()
  # not sure if we have to check the SOURCE_DIR for definitions
  get_target_property(_directory ${_name} SOURCE_DIR)
  get_property(
    _dir_defs
    DIRECTORY ${_directory}
    PROPERTY COMPILE_DEFINITIONS)
  if(_dir_defs)
    list(APPEND _target_defs ${_dir_defs})
  endif()
  # We could iterate trough all the sources to get their compile definitions,
  # but unfortunatelly COMPILE_DEFINITIONS are not visible outside the directory
  # where the target was added. Therefore the do not use this property.
  # get_target_property(_sources ${_name} SOURCES) foreach(_src ${_sources})
  # set(_src_loc ${_directory}/${_src}) get_source_file_property(_cpp_defs
  # "${_src_loc}" COMPILE_DEFINITIONS) if(_cpp_defs) list(APPEND _target_defs
  # "${_cpp_defs}") endif() endforeach()
  if(_target_defs)
    list(REMOVE_DUPLICATES _target_defs)
  endif()
  set(_defflags)
  foreach(_d ${_target_defs})
    set(_defflags ${_defflags} -D${_d})
  endforeach()
  set(${_output_name}
      "${_defflags}"
      PARENT_SCOPE)
endfunction()

# ==============================================================================
# FUNCTION: get_preprocessed_filename
# ==============================================================================
# Generate a name for the preprocessed source file _src_path is the the path of
# the source file with relative or absolute path _output_name is the name of the
# variable to store the results There is an optional suffix argument:
# get_preprocessed_filename(${source} _preproc_name "suffix")
function(get_preprocessed_filename _src_path _output_name)
  get_source_file_property(_loc ${_src_path} LOCATION) # I need the full path
  get_filename_component(_e "${_loc}" EXT)
  get_filename_component(_n "${_loc}" NAME_WE)
  get_filename_component(_dir "${_loc}" DIRECTORY)
  string(REGEX REPLACE "^${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" _dir
                       ${_dir})
  string(REGEX REPLACE "F" "f" _e "${_e}")
  if(${ARGC} GREATER 2)
    # add optional suffix
    set(${_output_name}
        "${_dir}/${_n}${ARGV2}${_e}"
        PARENT_SCOPE)
  else()
    set(${_output_name}
        "${_dir}/${_n}${_e}"
        PARENT_SCOPE)
  endif()
endfunction()
