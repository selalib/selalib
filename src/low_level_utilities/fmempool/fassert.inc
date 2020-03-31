!----------------------------------------------------------------------------
! C-like assert implementation for Fortran. To be used together with
! 'fassert.F90'.
!
! This file is part of the 'fmempool' package.  See README.md,
! LICENSE.txt, and the documentation for details.
!----------------------------------------------------------------------------

! https://gcc.gnu.org/onlinedocs/cpp/Traditional-macros.html#Traditional-macros
#if defined( __GFORTRAN__ )
#define CAT(x) "x"
#else
#define CAT(x) #x
#endif

#define ASSERT(x) \
if (.not. (x)) then; \
call assertion( CAT(x), __FILE__, __LINE__ ); \
end if;
