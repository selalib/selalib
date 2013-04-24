ADD_LIBRARY(sll_visu_pic STATIC sll_visu_pic.F90)

ADD_DEPENDENCIES(sll_visu_pic  sll_assertion 
                               sll_constants
                               sll_memory
                               sll_working_precision
                               sll_utilities
                               sll_file_io)

ADD_EXECUTABLE(test_visu_pic unit_test.F90 biot_savart.F90)

TARGET_LINK_LIBRARIES(test_visu_pic sll_visu_pic
                                    sll_file_io
                                    sll_constants
                                    sll_memory
                                    sll_assertion
                                    sll_working_precision
                                    sll_utilities)
