ADD_LIBRARY(sll_lin_trans_pic_2d STATIC sll_lin_trans_pic_2d.F90
                                        flow_base.F90
                                        leap_frog_1st_flow_vp_2d.F90
                                        leap_frog_2nd_flow_vp_2d.F90
                                        )

ADD_DEPENDENCIES(sll_lin_trans_pic_2d   sll_memory                                      
                                        sll_assertion 
                                        sll_working_precision
                                        )


#ADD_LIBRARY(sll_flow_base STATIC        flow_base.F90
#                                        leap_frog_1st_flow_vp_2d.F90
#                                        leap_frog_2nd_flow_vp_2d.F90
#                                        )
#


