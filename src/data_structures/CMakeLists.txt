ADD_SUBDIRECTORY(boundary_condition_descriptors)
ADD_SUBDIRECTORY(descriptors)
ADD_SUBDIRECTORY(distribution_function)
ADD_SUBDIRECTORY(fields)

IF(F2008)
  ADD_SUBDIRECTORY( vector_space )
ENDIF(F2008)
