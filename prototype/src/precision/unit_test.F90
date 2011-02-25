program working_precision_tester
#include "sll_working_precision.h"
  implicit none

  sll_int32   :: i =  z'7fffffff' !2147483647, largest 32-bit int
  sll_int64   :: i2
  sll_real32  :: f = z'7f7fffff'     ! largest IEEE-754 single prec float
  sll_real32  :: small = z'7f000001' ! smallest num with same binary exp as f
  sll_real32  :: neglected = z'7effffff'
  sll_real64  :: f2 = z'7fefffffffffffff' ! largest IEEE-754 double float
  sll_real64  :: small2 = z'7feffffffffffff1' 

  print *, '*************************************'
  print *, 'Tester for the working precision module'
  print *,
  print *, '*************************************'
  print *, 'Test the 32-bit int'
  write (*, '("  value of the largest representable i32 is ", i12)') i
  write (*, '("  KIND number for our i32 integer is", I2)') kind(i)
  write (*, '("  value of the largest i32 + 1 is ", i20)') i+1
  print *,
  print *, '*************************************'
  print *, 'Test the 32-bit float'
  write (*, '("  value of the largest representable f32 is ", es20.10)') f
  write (*, '("  KIND number for f32 is", I2)') kind(f)
  write (*, '("  smallest number with the same exponent as largest f32: ", es20.10)') small
  write (*, '("  value of the largest f32 + small is ", es20.10)') f+small
  print *, ' ...If the last value is +Infinity, our f32 got overflowed, which was the intent.'
  write (*, '("  largest number that will be neglected when added to the largest f32: ", es20.10)') neglected
  write (*, '("  value of the largest f32 + neglected is ", es20.10)') f+neglected
  print *,
  print *, '*************************************'
  write (*, '("  KIND number for our i64 integer is", I2)') kind(i2)
  print *, ' compare the last number with kind number for i32'
  write (*, '("  value of the largest representable f64 is ", es20.10)') f2
  write (*, '("  the value of the smallest f64 with same binary exp ", es20.10)') small2
  write (*, '("  value of the largest f32 + small is ", es20.10)') f2+small2
  print *,


  write ( *, '("KIND number for default integer is", I2)') kind ( 0 )
  write ( *, '("KIND number for i_8 integer is", I2)') kind(0_8)
  write ( *, '("KIND number for default single precision is", I2)') kind(0.0)
  write ( *, '("KIND number for default double precision is", I2)') kind(0.0D0)
  write ( *, '("KIND number for 2147483647 (2^32-1) is", I2)') kind(2147483647)
  write ( *, '("KIND number for 2147483648 (2^32) is", I2)')   kind(1_sll_i64)

contains

end program working_precision_tester
