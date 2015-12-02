ADD_SUBDIRECTORY(boundary_condition_descriptors)
ADD_SUBDIRECTORY(descriptors)
ADD_SUBDIRECTORY(distribution_function)
ADD_SUBDIRECTORY(fields)

IF(FULL_FORTRAN2003)
  ADD_SUBDIRECTORY( vector_space )
ENDIF(FULL_FORTRAN2003)
