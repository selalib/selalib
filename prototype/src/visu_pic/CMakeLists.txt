ADD_LIBRARY(sll_visu_pic STATIC sll_visu_pic.F90)

ADD_DEPENDENCIES(sll_visu_pic  sll_assertion 
                               sll_working_precision
                               sll_low_level_file_io)
