MESSAGE(STATUS "IBM system using xlf")

######################################################################
#AIX_BIT_FLAGS is the option to switch between 32/64
#If BITS=64, compiler and ar options have to have -q64 or -X64
######################################################################
SET(AIX_BIT_FLAGS "")
SET(F77OPTFLAGS  -O3 -qstrict -q64)
SET(AR_OPTIONS "-X64")
SET(F77 xlf)
SET(F77OPTFLAGS  -O3 -qstrict)
SET(FORTRAN_LIBS " -lxlf90 -lxlf")
SET(F77FLAGS ${F77OPTFLAGS})
