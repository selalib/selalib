ADD_LIBRARY(sll_coordinate_transformations STATIC
            sll_m_coordinate_transformation_2d_base.F90
            sll_m_coordinate_transformations_2d.F90
            sll_m_common_coordinate_transformations.F90)

TARGET_LINK_LIBRARIES(sll_coordinate_transformations 
                      sll_interpolators sll_file_io sll_meshes)
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
