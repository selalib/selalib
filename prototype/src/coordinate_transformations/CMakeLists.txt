FIND_PROGRAM(PYTHON_EXECUTABLE NAMES python3 python3.3 DOC "python")

FILE(GLOB TXT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.txt")
FILE(COPY ${TXT_FILES} DESTINATION ${CMAKE_BINARY_DIR})

IF(PYTHON_EXECUTABLE)

   SET(TRANSLATE_PYTHON_SCRIPT 
      "${CMAKE_CURRENT_SOURCE_DIR}/../python_scripts/translate_multipatch_info.py")

   SET(TRANSLATE_NML_FILE "${CMAKE_CURRENT_SOURCE_DIR}/circle_5mp_info.txt")

   SET(NURBS_PATCH_PYTHON_SCRIPT 
         "${CMAKE_CURRENT_SOURCE_DIR}/../python_scripts/nurbs_patch_txt_to_nml.py")

   FILE(GLOB CIRCLE_TXT_FILES "${CMAKE_CURRENT_SOURCE_DIR}/circle_5mp_patch*.txt")

   FOREACH(_TXT_CIRCLE_FILE ${CIRCLE_TXT_FILES})
      EXECUTE_PROCESS(
         COMMAND ${PYTHON_EXECUTABLE} ${NURBS_PATCH_PYTHON_SCRIPT} ${_TXT_CIRCLE_FILE}
      )
   ENDFOREACH()

   EXECUTE_PROCESS(
      COMMAND ${PYTHON_EXECUTABLE}  ${TRANSLATE_PYTHON_SCRIPT} ${TRANSLATE_NML_FILE}
      WORKING_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_FILE "circle_patch_info.out")

ELSE(PYTHON_EXECUTABLE)

   MESSAGE(AUTHOR_WARNING "Multipatch test needs python3")
   SET(PYTHON3_FOUND OFF)

ENDIF(PYTHON_EXECUTABLE)


#MESSAGE(STATUS "Coordinate transformation")
#MESSAGE(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
#MESSAGE(STATUS ${CMAKE_BINARY_DIR})
#MESSAGE(STATUS ${TXT_FILES})

FILE(GLOB NML_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.nml")
FILE(COPY ${NML_FILES} DESTINATION ${CMAKE_BINARY_DIR})

ADD_LIBRARY(sll_coordinate_transformations STATIC
            sll_coordinate_transformation_2d_base.F90
            sll_coordinate_transformations_2d.F90
            sll_coordinate_transformations_2d_nurbs.F90
            sll_common_coordinate_transformations.F90)

TARGET_LINK_LIBRARIES(sll_coordinate_transformations sll_interpolators
	                                               sll_file_io
	                                               sll_logical_meshes)

IF(HDF5_IS_PARALLEL)
   SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)

ADD_EXECUTABLE(test_coordinate_transformations_2d unit_test_2d.F90)

TARGET_LINK_LIBRARIES(test_coordinate_transformations_2d 
                         sll_coordinate_transformations)

