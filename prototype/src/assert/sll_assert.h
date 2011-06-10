#ifndef _assert_h_
#define _assert_h_

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


    ! Note the semicolon when the SLL_ASSERT() macro gets expanded. We need
    ! that this be a part of the macro since we also want to be able to
    ! use this macro call within other macros. If the expansion yields 
    ! nothing (second case) then we don't want dangling semicolons...
#ifdef DEBUG
# define SLL_ASSERT(x) if ( .not. (x) ) \
      call sll_assert( STRNG(x), __FILE__, __LINE__ );
#else
# define SLL_ASSERT(x) 
#endif

use sll_assertion


#endif
