#ifndef _SLL_REMAP_H_
#define _SLL_REMAP_H_

#include "misc_utils.h"

use remapper

! The intent of this macro is to get rid of the extremely ugly fact that
! the call to new_remap_plan_XD() requires also a call to the INT32_SIZEOF()
! macro. The solution is to pass either the array or an element of the 
! array for which the plan is needed... but in doing so, we would make this
! call type-dependent. Hence the macro is necessary in this case.

! Unfortunately, the following macro cannot be made typeless because it
! implicitly depends on the dimension of the passed array.
#define NEW_REMAP_PLAN_2D( arg1, arg2, array2d ) \
   new_remap_plan_2d( arg1, arg2, INT32_SIZEOF(array2d(1,1)) )

#define NEW_REMAP_PLAN_3D( arg1, arg2, array3D ) \
   new_remap_plan_3D( arg1, arg2, INT32_SIZEOF(array3D(1,1,1)) )

#define NEW_REMAP_PLAN_4D( arg1, arg2, array4D ) \
   new_remap_plan_4D( arg1, arg2, INT32_SIZEOF(array4D(1,1,1,1)) )

#define NEW_REMAP_PLAN_6D( arg1, arg2, array6D ) \
  new_remap_plan_6D( arg1, arg2, INT32_SIZEOF(array6D(1,1,1,1,1,1)) )


#endif
