#ifndef _assert_h_
#define _assert_h_

#include "misc_utils.h" 

    ! Note the semicolon when the SLL_ASSERT() macro gets expanded. We need
    ! that this be a part of the macro since we also want to be able to
    ! use this macro call within other macros. If the expansion yields 
    ! nothing (second case) then we don't want dangling semicolons...
#ifdef DEBUG
# define SLL_ASSERT(x) if ( .not. (x) ) then;          \
      call sll_assert( STRNG(x), __FILE__, __LINE__ ); \
   end if;
#else
# define SLL_ASSERT(x) 
#endif

use sll_assertion


#endif
