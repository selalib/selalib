SET(PIC2D_SRC zone.F90 
              quietstart.F90 
              particules.F90
              diagno.F90 )

ADD_LIBRARY(pic_vm_2d2v_cart STATIC ${PIC2D_SRC})

TARGET_LINK_LIBRARIES(pic_vm_2d2v_cart sll_pic_visu nufft)

ADD_EXECUTABLE(sim_pic_vm_2d2v_cart main.F90)
TARGET_LINK_LIBRARIES(sim_pic_vm_2d2v_cart pic_vm_2d2v_cart sll_poisson_solvers) 

#SET(ARGS "${CMAKE_CURRENT_SOURCE_DIR}/input/data.nml")
#ADD_TEST(NAME test_pic2d COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_pic2d ${ARGS})
#ET_TESTS_PROPERTIES(test_pic2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
