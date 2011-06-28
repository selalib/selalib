#ifndef _sll_splines_h_
#define _sll_splines_h_

#include "sll_assert.h"

  use sll_splines
  
  ! ***************************************************************************
  ! This header file exposes the available functionality for the spline-
  ! related types. As a user of this file one only needs to know that there
  ! is a type called 'sll_spline_1D' that can be:
  ! - declared
  ! - allocated (with new_spline_1D())
  ! - initialized
  ! - queried for its atributes through access functions (in this case, macros).
  !
  ! ***************************************************************************

#define PERIODIC_SPLINE   0_i32
#define NATURAL_SPLINE    1_i32

#define GET_SPLINE_DELTA(obj)    obj%delta
#define GET_SPLINE_XMIN(obj)     obj%xmin
#define GET_SPLINE_XMAX(obj)     obj%xmax





#endif
