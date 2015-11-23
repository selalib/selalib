ADD_LIBRARY( pppack STATIC 
  banfac.f90 
  banslv.f90  
  bsplvb.f90
)

ADD_SUBDIRECTORY(testing)
