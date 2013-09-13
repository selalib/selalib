ADD_LIBRARY(sll_remap STATIC sll_remap.F90)

TARGET_LINK_LIBRARIES(sll_remap sll_collective 
                                sll_electric_field_accumulators)  # this has to disappear


ADD_EXECUTABLE(test_remap_2d unit_test_2d.F90)
TARGET_LINK_LIBRARIES(test_remap_2d sll_remap)

ADD_EXECUTABLE(test_remap_3d unit_test_3d.F90)
TARGET_LINK_LIBRARIES(test_remap_3d sll_remap)

ADD_EXECUTABLE(test_remap_4d unit_test_4d.F90)
TARGET_LINK_LIBRARIES(test_remap_4d sll_remap)

ADD_EXECUTABLE(test_remap_6d unit_test_6d.F90)
TARGET_LINK_LIBRARIES(test_remap_6d sll_remap)

IF(HDF5_IS_PARALLEL AND HDF5_PARALLEL_ENABLED)

   ADD_EXECUTABLE(test_layout_output unit_test_output.F90)
   TARGET_LINK_LIBRARIES(test_layout_output sll_remap sll_file_io_parallel)

ENDIF()



