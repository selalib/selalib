# This step is only necessary when we want the .txt files to be available in the
# build directory for use by the unit tests. Geometric files aren't needed, only
# their .nml translations. Hence this is commented out. FILE(GLOB TXT_FILES
# "${CMAKE_CURRENT_SOURCE_DIR}/*.txt") FILE(COPY ${TXT_FILES} DESTINATION
# ${CMAKE_BINARY_DIR})

if(PYTHON3_FOUND)

  # choose python script to translate individual patch information.
  set(PATCH_SCRIPT "${CMAKE_SOURCE_DIR}/python/nurbs_patch_txt_to_nml.py")

  # Search through subdirectories to do the translations of the patch and
  # connectivity information.
  file(GLOB_RECURSE PATCH_FILES FOLLOW_SYMLINKS
       "${CMAKE_CURRENT_SOURCE_DIR}/*_patch*.txt")
  # MESSAGE(${PATCH_FILES})

  # make a call to nurbs_patch_txt_to_nml.py something_patchi.txt in each patch
  # file found in the subdirectories. Since each file found by CMake is
  # prepended with its path, the output also should be written in the same
  # directory.
  foreach(PATCH ${PATCH_FILES})
    # MESSAGE(STATUS "converting file:" ${PATCH})
    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${PATCH_SCRIPT} ${PATCH})
    unset(PATCH)
  endforeach()

  set(MP_INFO_SCRIPT "${CMAKE_SOURCE_DIR}/python/translate_multipatch_info.py")

  file(GLOB_RECURSE INFO_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*_info.txt")

  foreach(PATCH ${INFO_FILES})
    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${MP_INFO_SCRIPT} ${PATCH})
    unset(PATCH)
  endforeach()

  # Algorithm to generate the files needed for a multipatch representation: 1.
  # Find if the file ID.txt is available 2. If ID.txt is available, find the
  # 'root-name' of the transformation by extracting it from whatever file that
  # contains the "_info.txt" ending. 3. Count how many patches are there in the
  # transformation. Points 2. and 3. are carried out by the Python script
  # itself. 4. Assemble the call and call the appropriate script with the
  # desired arguments.

  set(CONNECTIVITY_SCRIPT
      "${CMAKE_SOURCE_DIR}/python/translate_connectivity_info.py")

  # we had read already the *_info.txt in INFO_FILES. FILE(GLOB_RECURSE
  # IEN_TXT_FILES FOLLOW_SYMLINKS "${CMAKE_CURRENT_SOURCE_DIR}/IEN_*.txt")
  # FILE(GLOB_RECURSE LM_TXT_FILES FOLLOW_SYMLINKS
  # "${CMAKE_CURRENT_SOURCE_DIR}/LM_*.txt")

  file(GLOB_RECURSE ID_TXT_FILES FOLLOW_SYMLINKS
       "${CMAKE_CURRENT_SOURCE_DIR}/ID.txt")

  # list(LENGTH IEN_TXT_FILES NUM_PATCHES) MESSAGE("${IEN_TXT_FILES}")
  # MESSAGE("Processing a multipatch with ${NUM_PATCHES} patches.")

  foreach(ID_FILE ${ID_TXT_FILES})
    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CONNECTIVITY_SCRIPT}
                            ${ID_FILE})
  endforeach()

else(PYTHON3_FOUND)
  message(STATUS "Multipatch test needs python3")
endif(PYTHON3_FOUND)

# .nml files usually need to be accessed by the unit tests, so copy them to the
# build directory. YG  FILE(GLOB_RECURSE NML_FILES
# "${CMAKE_CURRENT_SOURCE_DIR}/*.nml") YG  FILE(COPY ${NML_FILES} DESTINATION
# ${CMAKE_BINARY_DIR})
