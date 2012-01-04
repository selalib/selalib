program unit_test
  use sll_fft
  use numeric_constants
  use sll_timer

#include "sll_working_precision.h"
  implicit none
  
  sll_int32 :: i
  sll_comp64, allocatable, dimension(:)  :: dat, in
  sll_comp32, allocatable, dimension(:)  :: c
  type(time_mark), pointer :: mark 
  sll_real64 :: time, time2
  sll_int32 :: s, n
  type(sll_fft_plan), pointer :: fft, fft2, fft3  
  mark => new_time_mark() 

  do s=16,20
    n=2**s
    print *,'n=',n
    allocate(dat(0:n-1))
    allocate(in(0:n-1))
    allocate(c(0:n-1))

    dat(0) = 1.d0
    dat(1:n-1) = 0.d0

    in(0) = 1.d0
    in(1:n-1) = 0.d0

    c(0) = 1.d0
    c(1:n-1) = 0.d0

    mark => reset_time_mark(mark)
    fft => sll_new_fft(n,FFT_COMPLEX)
    !time2 = time_elapsed_since(mark)
    !print *, 'SLL_FFT : initializing time : ',time2
    !mark => reset_time_mark(mark)
    fft => sll_apply_fft(fft,dat,FFT_INVERSE)
    time2 = time_elapsed_since(mark)
    print *, 'SLL_FFT : fft time : ',time2

    mark => start_time_mark(mark)
    fft2 => sll_new_fft(n,FFT_COMPLEX,FFTW_MOD)
    !time = time_elapsed_since(mark)
    !print *, 'FFT_PACK : initializing time : ',time
    !mark => reset_time_mark(mark)
    fft2 => sll_apply_fft(fft2,in,FFT_INVERSE)
    time = time_elapsed_since(mark)
    print *, 'FFTW : fft time : ',time

    mark => reset_time_mark(mark)
    fft3 => sll_new_fft(n,FFT_COMPLEX,FFTPACK_MOD)
    !time2 = time_elapsed_since(mark)
    !print *, 'SLL_FFT : initializing time : ',time2
    !mark => reset_time_mark(mark)
    fft3 => sll_apply_fft(fft3,c,FFT_INVERSE)
    time2 = time_elapsed_since(mark)
    print *, 'FFTPACK : fft time : ',time2


    !open(4,file='dat.txt')
    !do i=0,n-1
    !  write(4,*) dat(i) , in(i), c(i)
    !enddo  
    !close(4)

    !print *, MAXVAL(ABS(2*in-dat-c))
    deallocate(dat)
    deallocate(in)
    deallocate(c)
    fft => sll_delete_fft(fft)  
    fft2 => sll_delete_fft(fft2)  
    fft3 => sll_delete_fft(fft3)  
  enddo
  
end program unit_test
