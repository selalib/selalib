#ifndef _SLL_WORKING_PRECISION_
#define _SLL_WORKING_PRECISION_

  ! For the definition of the kinds, refer to the file 
  ! sll_working_precision.F90
  !
  ! We provide a few aliases for convenience. Whenever multiple precision
  ! representations or computations are needed, these aliases provide a
  ! centralized place to define the numeric type. As their name suggests,
  ! the intent of something like sll_int64 is to guarantee a 64-bit integer,
  ! likewise for the other aliases.
  !
  ! We also provide the unqualified sll_int and sll_real which can be used
  ! anywhere, while providing a single location for changing the size of the
  ! representation library-wide.
  ! 
  ! To qualify numerical constants, like 2011_i32 we offer _i32, _i64,
  ! _f32, and _f64 for the moment.

  ! These are left without qualification and meant as the 'user' type
#define sll_int    integer(kind=i32)
#define sll_real   real(kind=f64)

  ! Selalib's standard aliases:

  ! Integer types:
#define sll_int32  integer(kind=i32)
#define sll_int64  integer(kind=i64)
  ! Floating point types:
#define sll_real32 real(kind=f32)
#define sll_real64 real(kind=f64)
  
  ! Complex types
#define sll_comp32 complex(kind=f32)
#define sll_comp64 complex(kind=f64)

use sll_working_precision

#endif

