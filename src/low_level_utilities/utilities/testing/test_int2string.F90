program test_int2_string
#include "sll_working_precision.h"
use sll_m_utilities
implicit none

sll_int32 :: istep, j
character(len=1) :: cstep1
character(len=2) :: cstep2
character(len=3) :: cstep3
character(len=4) :: cstep4
character(len=5) :: cstep5

do j = 1, 10
  istep = 2**j
  call sll_s_int2string(istep, cstep1)
  write(*,*) cstep1
  istep = 99
  call sll_s_int2string(istep, cstep2)
  write(*,*) cstep2
  istep = 999
  call sll_s_int2string(istep, cstep3)
  write(*,*) cstep3
  istep = 9999
  call sll_s_int2string(istep, cstep4)
  write(*,*) cstep4
  istep = 99999
  call sll_s_int2string(istep, cstep5)
  write(*,*) cstep5
end do


end program test_int2_string
