#ifndef _sll_interpolators_h_
#define _sll_interpolators_h_

#ifndef STDF95
use sll_module_interpolators_1d_base
#endif
use sll_cubic_spline_interpolator_1d
use sll_quintic_spline_interpolator_1d
use sll_odd_degree_spline_interpolator_1d

#ifndef STDF95
use sll_module_interpolators_2d_base
#endif
use sll_cubic_spline_interpolator_2d

#include "sll_splines.h"

#endif
 
