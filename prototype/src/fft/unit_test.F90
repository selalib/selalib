program test
  use sll_fft
  !use sll_timer
#include "sll_working_precision.h"
  implicit none

  sll_int32, parameter :: n = 2**3
  sll_int32, parameter :: m = 2**2
  type(sll_fft_plan), pointer :: plan
  type(sll_fft_plan_2d), pointer :: plan2d
  sll_comp64, dimension(n) :: data_comp
  sll_comp64, dimension(n,m) :: data_comp_2d
  sll_real64, dimension(n) :: data_real
  !sll_comp64, dimension(4,2,2) :: data2
  sll_int32 :: i,j

  print *,'-------------------------------------------------'
  print * ,'COMPLEX TO COMPLEX'

  data_comp(1) = complex(1_f64,0_f64)
  data_comp(2:) = complex(0_f64,0_f64)
  
  print *,''
  print *,'input :'
  do i=1,n
    print *, data_comp(i)
  enddo
  print *,'perform forward normalize' 
  plan => new_plan_c2c_1d(n,data_comp,data_comp,FFT_FORWARD,FFT_NORMALIZE_FORWARD)
  call apply_fft_c2c_1d(plan,data_comp,data_comp)
  call delete(plan)
  do i=1,n
    print *, data_comp(i)
  enddo

  print *,'perform inverse' 
  plan => new_plan_c2c_1d(n,data_comp,data_comp,FFT_INVERSE)
  call apply_fft_c2c_1d(plan,data_comp,data_comp)
  call delete(plan)
  do i=1,n
    print *, data_comp(i)
  enddo

  print *,'perform forward' 
  plan => new_plan_c2c_1d(n,data_comp,data_comp,FFT_FORWARD)
  call apply_fft_c2c_1d(plan,data_comp,data_comp)
  call delete(plan)
  do i=1,n
    print *, data_comp(i)
  enddo

  print *,'perform inverse normalize' 
  plan => new_plan_c2c_1d(n,data_comp,data_comp,FFT_INVERSE,FFT_NORMALIZE_INVERSE)
  call apply_fft_c2c_1d(plan,data_comp,data_comp)
  call delete(plan)
  do i=1,n
    print *, data_comp(i)
  enddo


  print * ,''
  print *,'-------------------------------------------------'
  print * ,''
  print * ,'REAL TO REAL'
  print * ,''

  data_real(1) = 1_f64
  data_real(2:) = 0_f64
 
  print *,''
  print *,'input :'
  do i=1,n
    print *, data_real(i)
  enddo

  print *,'perform r2r forward'
  plan => new_plan_r2r_1d(n,data_real,data_real,FFT_FORWARD,FFT_NORMALIZE_FORWARD)
  call apply_plan_r2r_1d(plan,data_real,data_real)
  call delete(plan)
  print *,'perform r2r inverse'
  plan => new_plan_r2r_1d(n,data_real,data_real,FFT_INVERSE)
  call apply_plan_r2r_1d(plan,data_real,data_real)
  call delete(plan)

  print *,''
  do i=1,n
    print *, data_real(i)
  enddo

  print *,'-------------------------------------------------'
  print * ,''
  print * ,'REAL TO COMPLEX'
  print * ,''
  data_real(1) = 1_f64
  data_real(2:) = 0_f64
  print *,'input :'
  do i=1,n
    print *, data_real(i)
  enddo

  plan => new_plan_r2c_1d(n,data_real,data_comp)
  call apply_plan_r2c_1d(plan,data_real,data_comp)
  call delete(plan)

  print * ,''
  do i=1,n/2+1
    print *, data_comp(i)
  enddo

  print * ,''
  print * ,'COMPLEX TO REAL'
  print * ,''
  plan => new_plan_c2r_1d(n,data_comp(1:n/2+1),data_real,FFT_NORMALIZE_INVERSE)
  call apply_plan_c2r_1d(plan,data_comp(1:n/2+1),data_real)
  call delete(plan)
  do i=1,n
    print *, data_real(i)
  enddo

  print * ,''
  print *,'-------------------------------------------------'
  print * ,'COMPLEX TO COMPLEX 2D'

  data_comp_2d(:,:) = complex(0_f64,0_f64)
  data_comp_2d(1,:) = complex(1_f64,0_f64)
  
  print *,''
  print *,'input :'
  do i=1,n
   do j=1,m
    print * , i,j,data_comp_2d(i,j)
   enddo
  enddo
  print *,'perform forward'
  plan2d => new_plan_c2c_2d(n,m,data_comp_2d,data_comp_2d,FFT_FORWARD)
  call apply_plan_c2c_2d(plan2d,data_comp_2d,data_comp_2d)
  call delete(plan2d)
  print *,'perform inverse'
  plan2d => new_plan_c2c_2d(n,m,data_comp_2d,data_comp_2d,FFT_INVERSE,FFT_NORMALIZE_INVERSE)
  call apply_plan_c2c_2d(plan2d,data_comp_2d,data_comp_2d)
  call delete(plan2d)
  do i=1,n
   do j=1,m
    print *, i,j,data_comp_2d(i,j)
   enddo
  enddo
  print *,'-------------------------------------------------'
  print * ,''

!  data2 = complex(1_f64,0_f64)
 
 
 
!  plan => new_plan_c2c_1d(4,2,2,in,in,FFT_FORWARD_X)
!  call apply_fft_c2c_1d(plan,in,in)
!  call delete(plan)
!  print *, in(:,1,1)
!
!
!  plan => new_plan_c2c_1d(4,in,in,FFT_INVERSE)
!  call apply_fft_c2c_1d(plan,in,in)
!  call delete(plan)
!
!  do i=1,4
!    print *, in(i)
!  enddo


!contains
 
end program test

