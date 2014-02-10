FIND_PROGRAM(PYTHON_EXECUTABLE NAMES python3 DOC "python")

FILE(GLOB TXT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.txt")
FILE(COPY ${TXT_FILES} DESTINATION ${CMAKE_BINARY_DIR})

IF(PYTHON_EXECUTABLE)

   SET(TRANSLATE_PYTHON_SCRIPT 
      "${CMAKE_CURRENT_SOURCE_DIR}/../python_scripts/translate_multipatch_info.py")

   SET(TRANSLATE_NML_FILE "${CMAKE_CURRENT_SOURCE_DIR}/identity_mp_info.txt")

   SET(NURBS_PATCH_PYTHON_SCRIPT 
         "${CMAKE_CURRENT_SOURCE_DIR}/../python_scripts/nurbs_patch_txt_to_nml.py")

   FILE(GLOB IDENTITY_TXT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/identity_mp_patch*.txt")

   FOREACH(_TXT_IDENTITY_FILE ${IDENTITY_TXT_FILES})
      EXECUTE_PROCESS(
         COMMAND ${PYTHON_EXECUTABLE} ${NURBS_PATCH_PYTHON_SCRIPT} ${_TXT_IDENTITY_FILE}
	 WORKING_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}
         OUTPUT_FILE "nurbs_patch.out"
      )
   ENDFOREACH()

   EXECUTE_PROCESS(
      COMMAND ${PYTHON_EXECUTABLE}  ${TRANSLATE_PYTHON_SCRIPT} ${TRANSLATE_NML_FILE}
      WORKING_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_FILE "identity_patch_info.out")

ELSE(PYTHON_EXECUTABLE)

   MESSAGE(SEND_WARN "Multipatch test needs python3")

ENDIF(PYTHON_EXECUTABLE)


#MESSAGE(STATUS "Coordinate transformation multipatchs")
#MESSAGE(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
#MESSAGE(STATUS ${CMAKE_BINARY_DIR})
#MESSAGE(STATUS ${TXT_FILES})

FILE(GLOB NML_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.nml")
FILE(COPY ${NML_FILES} DESTINATION ${CMAKE_BINARY_DIR})


IF(NOT STDF95)

   ADD_LIBRARY(sll_coordinate_transformation_multipatch
     STATIC sll_coordinate_transformation_multipatch.F90)

   TARGET_LINK_LIBRARIES(sll_coordinate_transformation_multipatch
     sll_coordinate_transformations
     sll_logical_meshes_multipatch
     sll_memory)

   ADD_EXECUTABLE(test_coordinate_transformation_multipatch_2d unit_test_2d.F90)
   TARGET_LINK_LIBRARIES(test_coordinate_transformation_multipatch_2d
     sll_coordinate_transformation_multipatch
     sll_file_io
     sll_memory)

ENDIF(NOT STDF95)
