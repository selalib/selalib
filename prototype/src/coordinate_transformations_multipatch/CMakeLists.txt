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
