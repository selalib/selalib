program working_precision_tester
#include "sll_working_precision.h"
  ! use, intrinsic :: ieee_arithmetic ! Why is this module not loading???
  ! because gfortran is not standard compliant... 
  implicit none
 
  ! **************************************************************************
  ! Basically all we want to know is if the integer values have representations
  ! in 32 and 64 bits and if the same holds for the reals. Therefore, if the
  ! real kind parameter actually represents the number of bytes used for the
  ! storage of those numbers (4 & 8), we are done. The problem is that 
  ! according to the standard, that may not be the case (so the approach isn't 
  ! fully portable); here we just gather information that might be useful in
  ! learning the floating point features of such systems, and this was also
  ! an excuse to see what Fortran allowed us to do in this domain, as this
  ! might be useful for other purposes. Unfortunately, gfortran does not
  ! support some interesting ieee functionality that is supposed to be in a
  ! standard-compliant implementation.
  !
  ! Presently this behaves more like a playground than an actual tester. This
  ! should improve a lot more.
  !
  ! **************************************************************************

  sll_int32   :: test_int1 = 0
  sll_int32   :: test_int2 = 0
  sll_int32   :: i =  z'7fffffff' !2147483647, largest 32-bit int
  sll_int32   :: counter = 1
  sll_int32   :: exponent = -127
  sll_int64   :: i2
  sll_real32  :: test_float = 0.0
  sll_real32  :: ref_float = 0.0
  sll_real32  :: stepf32
  sll_real32  :: f = z'7f7fffff'     ! largest IEEE-754 single prec float
  sll_real32  :: small = z'7f000001' ! smallest num with same binary exp as f
  sll_real32  :: neglected =  z'72ffffff' ! 25-power of 2 orders less
  sll_real64  :: f2 = z'7fefffffffffffff' ! largest IEEE-754 double float
  sll_real64  :: small2 = z'7feffffffffffff1' 

  print *, '*************************************'
  print *, 'Tester for the working precision module'
  print *,
  print *, '*************************************'
  print *, 'Test the 32-bit int'
  ! test based on the change of sign that the integer suffers when there is
  ! an overflow
  infinite:  do 
     test_int1 = test_int1 + 1
     if (test_int1 .le. 0) then
        write (*, '(a, i12)') '  By measurement, the largest i32 is ', test_int2
        exit infinite
     else
        test_int2 = test_int2 + 1
     end if
  end do infinite
  write (*, '("  In theory, the value of the largest i32 is ", i12)') i
  write (*, '("  KIND number for our i32 integer is", I2)') kind(i)
  write (*, '("  value of the largest (theoretical) i32 + 1 is ", i20)') i+1
  print *,
  print *, '*************************************'
  print *, 'Test the 32-bit float'
  stepf32 = z'00000001' ! smallest positive normalized number
  infinite_real: do
     test_float = test_float + stepf32
     if ( test_float .gt. ref_float ) then
        counter = counter + 1
        ref_float = ref_float + stepf32
     else if (ref_float .eq. test_float) then ! yes, we compare the actual bits
        write (*, '(a, i12, a, i12)') '  Went through ', counter, &
             ' numbers for exponent ', exponent  
        stepf32 = stepf32*2.0   ! need the next bigger step
        exponent = exponent + 1 ! and move to the next exponent
        write (*, '(a, es20.12)') '  New step is ', stepf32
        counter = 1 ! reset the counter
        write (*, '(a, z15)') ' Largest number thus far: ', test_float
        if ( test_float .eq. transfer(z'7F800000', test_float) ) then
           exit infinite_real
        end if
     end if
  end do infinite_real
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
  write ( *, '("KIND number for 2147483648 (2^32) is", I2)')   kind(1_i64)

  write ( *, '(a, i4)') 'selected_int_kind(0): ', selected_int_kind(0)
  write ( *, '(a, i4)') 'size of sll_int64 is: ',  SLL_SIZEOF(i2)
print *, bit_size(counter)
contains

end program working_precision_tester
