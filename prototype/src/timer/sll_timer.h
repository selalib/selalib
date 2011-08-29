#ifndef _sll_time_h_
#define _sll_time_h_

! 'convenience' header file intended to hide the C-language dependence
! of the timing facility.

use sll_timer
use iso_c_binding

#define time_mark type(c_ptr)


#endif
