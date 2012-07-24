#ifndef _SLL_MISC_UTILS_H_
#define _SLL_MISC_UTILS_H_
  ! Some useful services that have not found a home in more specific
  ! modules.
  
  ! ************************************************************************
  ! Unfortunately, fpp does not recognize the cpp operators # and ##. Some
  ! compilers do, so we need this ugly workaround... and this only fixes
  ! the lack of #. ## is still a dream.
  !
  ! ************************************************************************

#if (defined ( GFORTRAN ) || defined ( G95 ) || defined(MPIF90))
# define STRNG(x) "x"
#else
# define STRNG(x) #x
#endif

  ! The following is useless in fpp apparently. The workaround of using
  ! double quotes does not allow to expand a macro for subsequent conversion
  ! to a string. We leave this here as a testament to what would have been
  ! nice to have.
#define XSTRNG( x ) STRNG( x )




! BYTE_SIZEOF() uses the byte_size, which is defined in misc_utils.F90. This
! macro returns the size of 'var' measured in bytes.
! INT32_SIZEOF() uses i32, which is defined in the basic numeric types:
! sll_working_precision.F90. The use of the integer as a yardstick is a
! more natural choice in some contexts.

#define BYTE_SIZEOF( var )  size(transfer(var, (/1_byte_size/) ))
#define INT32_SIZEOF( var ) size(transfer(var, (/1_i32/)))

use sll_misc_utils


#endif
