#ifndef _SLL_MISC_UTILS_H_
#define _SLL_MISC_UTILS_H_
!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

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




! BYTE_SIZEOF() uses the byte_size, which is defined in sll_utilities.F90. This
! macro returns the size of 'var' measured in bytes.
! INT32_SIZEOF() uses i32, which is defined in the basic numeric types:
! sll_working_precision.F90. The use of the integer as a yardstick is a
! more natural choice in some contexts.

#define BYTE_SIZEOF( var )  size(transfer(var, (/1_byte_size/) ))
#define INT32_SIZEOF( var ) size(transfer(var, (/1_i32/)))

#define SWAP(A,B) A = A + B; B = A - B; A = A - B

use sll_utilities
use sll_tridiagonal

#endif
