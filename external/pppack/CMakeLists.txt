ADD_LIBRARY( pppack STATIC 
  banfac.f90 
  banslv.f90  
  bsplvb.f90
  woodbury.F90
)

ADD_SUBDIRECTORY(testing)
