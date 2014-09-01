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

  sll_int32   :: test_int1
  sll_int32   :: test_int2
  sll_int32   :: istep
  sll_int32   :: i
  sll_int32   :: counter
  sll_int32   :: exponent
  sll_int64   :: i2
  sll_real32  :: test_float
  sll_real32  :: ref_float
  sll_real32  :: stepf32
  sll_real32  :: f
  sll_real32  :: small
  sll_real32  :: neglected
  sll_real64  :: f2
  sll_real64  :: small2

  test_int1  = 0
  test_int2  = 0
  i          = transfer(z'7fffffff',i)  !2147483647, largest 32-bit int
  counter    = 1
  exponent   = -127
  test_float = 0.0
  ref_float  = 0.0
  f          = transfer(z'7f7fffff',f)
  small      = transfer(z'7f000001',small) ! smallest num with same binary exp as f
  neglected  = transfer(z'72ffffff',neglected) ! 25-power of 2 orders less
  f2         = transfer(z'7fefffffffffffff',f2) ! largest IEEE-754 double float
  small2     = transfer(z'7feffffffffffff1',small2) 



  print *, '*************************************'
  print *, 'Tester for the working precision module'
  print *, ' '
  print *, '*************************************'
  print *,'#i32=',i32
  print *,'#i64=',i64
  print *,'#f32=',f32
  print *,'#f64=',f64
  print *, 'Test the 32-bit int'
  
  ! new test not so slow
  ! test based on the change of sign that the integer suffers when there is
  ! an overflow
  test_int1 = 1
  do istep=1,32
    test_int2 = test_int2+test_int1
    if(test_int2 .le. 0) then
      print *,'#Problem 2**',istep,'-1 should be positive'      
      stop
    endif
    test_int1 = 2*test_int1
    if(test_int1 .le. 0)then
       print *,'#2**',istep,'-1=',test_int2
       write (*, '(a, i12)') '#  By measurement, the largest i32 is ', test_int2
       exit
    else    
      print *,'#2**',istep,'-1=',test_int2,'2**',istep,'=',test_int1
    endif   
  enddo
  
  ! old test which may be too slow
  ! test based on the change of sign that the integer suffers when there is
  ! an overflow
!  test_int1 = 0
!  test_int2 = 0
!  infinite:  do 
!     test_int1 = test_int1 + 1
!     if (test_int1 .le. 0) then
!        write (*, '(a, i12)') '  By measurement, the largest i32 is ', test_int2
!        exit infinite
!     else
!        test_int2 = test_int2 + 1
!     end if
!  end do infinite
  write (*, '("  In theory, the value of the largest i32 is ", i12)') i
  write (*, '("  KIND number for our i32 integer is", I2)') kind(i)
  write (*, '("  value of the largest (theoretical) i32 + 1 is ", i20)') i+1
  print *, ' '
  print *, '*************************************'
  print *, 'Test the 32-bit float'
  stepf32 = transfer(z'00000001',stepf32) ! smallest positive normalized number
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
        write (*, '(a, es20.12)') ' Largest number thus far: ', test_float
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
  print *, ' '
  print *, '*************************************'
  write (*, '("  KIND number for our i64 integer is", I2)') kind(i2)
  print *, ' compare the last number with kind number for i32'
  write (*, '("  value of the largest representable f64 is ", es20.10)') f2
  write (*, '("  the value of the smallest f64 with same binary exp ", es20.10)') small2
  write (*, '("  value of the largest f32 + small is ", es20.10)') f2+small2
  print *, ' '


  write ( *, '("KIND number for default integer is", I2)') kind ( 0 )
  write ( *, '("KIND number for i_8 integer is", I2)') kind(0_8)
  write ( *, '("KIND number for default single precision is", I2)') kind(0.0)
  write ( *, '("KIND number for default double precision is", I2)') kind(0.0D0)
  write ( *, '("KIND number for 2147483647 (2^32-1) is", I2)') kind(2147483647)
  write ( *, '("KIND number for 2147483648 (2^32) is", I2)')   kind(1_i64)

  write ( *, '(a, i4)') 'selected_int_kind(0): ', selected_int_kind(0)
print *, bit_size(counter)

end program working_precision_tester
