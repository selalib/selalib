SET(PIC2D_SRC zone.F90 
              initialisation.F90 
              quietstart.F90 
              particules.F90
	        maxwell.F90 
              diagno.F90 )

ADD_LIBRARY(pic2d STATIC ${PIC2D_SRC})

TARGET_LINK_LIBRARIES(pic2d sll_pic_visu nufft)

ADD_EXECUTABLE(test_pic2d main.F90)
TARGET_LINK_LIBRARIES(test_pic2d pic2d sll_poisson_solvers) 

SET(ARGS "${CMAKE_CURRENT_SOURCE_DIR}/input/data.nml")
ADD_TEST(NAME test_pic2d COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_pic2d ${ARGS})
SET_TESTS_PROPERTIES(test_pic2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
