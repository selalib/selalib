#ifndef _SLL_REMAP_H_
#define _SLL_REMAP_H_

#include "misc_utils.h"

use remapper

! The intent of this macro is to get rid of the extremely ugly fact that
! the call to new_remap_plan_XD() requires also a call to the INT32_SIZEOF()
! macro. The solution is to pass either the array or an element of the 
! array for which the plan is needed... but in doing so, we would make this
! call type-dependent. Hence the macro is necessary in this case.

#define NEW_REMAPPER_PLAN_3D( arg1, arg2, array3D ) \
   new_remap_plan_3D( arg1, arg2, INT32_SIZEOF(array3D(1,1,1)) )


#endif
