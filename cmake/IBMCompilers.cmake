message(STATUS "IBM system using xlf")

# ##############################################################################
# AIX_BIT_FLAGS is the option to switch between 32/64 If BITS=64, compiler and
# ar options have to have -q64 or -X64
# ##############################################################################
set(AIX_BIT_FLAGS "")
set(F77OPTFLAGS -O3 -qstrict -q64)
set(AR_OPTIONS "-X64")
set(F77 xlf)
set(F77OPTFLAGS -O3 -qstrict)
set(FORTRAN_LIBS " -lxlf90 -lxlf")
set(F77FLAGS ${F77OPTFLAGS})
