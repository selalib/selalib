#ifndef _assert_h_
#define _assert_h_

  ! ************************************************************************
  ! Unfortunately, fpp does not recognize the cpp directives # and ##. Some
  ! compilers do, so we need this ugly workaround... and this only fixes
  ! the lack of #. ## is still a dream.
  !
  ! ************************************************************************

#if (defined ( GFORTRAN ) || defined ( G95 ) || defined( PGI )
# define STRNG(x) "x"
#else
# define STRNG(x) #x
#endif
#define XSTRING(x) STRNG(x)
#ifndef NDEBUG
# define SLL_ASSERT(x, msg) if (.not. (x) ) \
  call sll_assert( STRNG(x)//" assertion failed in file "//__FILE__//", line: "//XSTRNG(__LINE__)//".")
#else
# define SLL_ASSERT(x, msg) 
#endif

use sll_assert


#endif
