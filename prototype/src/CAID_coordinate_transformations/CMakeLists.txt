FIND_PROGRAM(PYTHON_EXECUTABLE NAMES python3 python3.3 DOC "python")

# This step is only necessary when we want the .txt files to be available
# in the build directory for use by the unit tests. Geometric files aren't
# needed, only their .nml translations. Hence this is commented out.
#FILE(GLOB TXT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.txt")
#FILE(COPY ${TXT_FILES} DESTINATION ${CMAKE_BINARY_DIR})

IF(PYTHON_EXECUTABLE)

  SET(PYTHON3_FOUND ON)
  
  # choose python script to translate individual patch information.
  SET(PATCH_SCRIPT 
    "${CMAKE_CURRENT_SOURCE_DIR}/../python_scripts/nurbs_patch_txt_to_nml.py")
  
  # Search through subdirectories to do the translations of the patch and
  # connectivity information.  
  FILE(GLOB_RECURSE PATCH_FILES FOLLOW_SYMLINKS 
    "${CMAKE_CURRENT_SOURCE_DIR}/*_patch*.txt")
#  MESSAGE(${PATCH_FILES})
   

   # make a call to nurbs_patch_txt_to_nml.py something_patchi.txt in each
   # patch file found in the subdirectories. Since each file found by CMake
   # is prepended with its path, the output also should be written in the
   # same directory.
   FOREACH(PATCH ${PATCH_FILES})
     # MESSAGE("converting file:" ${PATCH})
     EXECUTE_PROCESS(
       COMMAND ${PYTHON_EXECUTABLE} ${PATCH_SCRIPT} ${PATCH}
     )
   ENDFOREACH()

   SET(MP_INFO_SCRIPT 
      "${CMAKE_CURRENT_SOURCE_DIR}/../python_scripts/translate_multipatch_info.py")
     
   FILE(GLOB_RECURSE INFO_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*_info.txt")

   FOREACH(PATCH ${INFO_FILES})	
      EXECUTE_PROCESS(
	COMMAND ${PYTHON_EXECUTABLE}  ${MP_INFO_SCRIPT} ${PATCH}
	)
   ENDFOREACH()
	  
   # Algorithm to generate the files needed for a multipatch representation:
   # 1. Find if the file ID.txt is available
   # 2. If ID.txt is available, find the 'root-name' of the transformation
   #    by extracting it from whatever file that contains the "_info.txt" 
   #    ending.
   # 3. Count how many patches are there in the transformation. Points 2. and
   #    3. are carried out by the Python script itself.
   # 4. Assemble the call and call the appropriate script with the desired
   #    arguments.

   SET(CONNECTIVITY_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/../python_scripts/translate_connectivity_info.py")

   # we had read already the *_info.txt in INFO_FILES.
#   FILE(GLOB_RECURSE IEN_TXT_FILES FOLLOW_SYMLINKS 
#     "${CMAKE_CURRENT_SOURCE_DIR}/IEN_*.txt")
#   FILE(GLOB_RECURSE LM_TXT_FILES FOLLOW_SYMLINKS 
#     "${CMAKE_CURRENT_SOURCE_DIR}/LM_*.txt")
   FILE(GLOB_RECURSE ID_TXT_FILES FOLLOW_SYMLINKS 
     "${CMAKE_CURRENT_SOURCE_DIR}/ID.txt")

#   list(LENGTH IEN_TXT_FILES NUM_PATCHES)
#   MESSAGE("${IEN_TXT_FILES}")
#   MESSAGE("Processing a multipatch with ${NUM_PATCHES} patches.")

   EXECUTE_PROCESS(
    	COMMAND ${CONNECTIVITY_SCRIPT} ${ID_TXT_FILES}
        )

ELSE(PYTHON_EXECUTABLE)
   MESSAGE(STATUS "Multipatch test needs python3")
   SET(PYTHON3_FOUND OFF)
ENDIF(PYTHON_EXECUTABLE)

# .nml files usually need to be accessed by the unit tests, so copy them to
# the build directory.
FILE(GLOB_RECURSE NML_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.nml")
FILE(COPY ${NML_FILES} DESTINATION ${CMAKE_BINARY_DIR})


