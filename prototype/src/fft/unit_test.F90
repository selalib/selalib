program unit_test
  use sll_fft
  use numeric_constants
  use sll_timer

#include "sll_working_precision.h"
  implicit none
  
  
  sll_int32 :: i
  sll_real64, allocatable, dimension(:)  :: copy_data_real, data_real_sll, data_real_w, data_real_pack
  type(time_mark), pointer :: mark 
  sll_real64 :: time, val, phase, acc
  sll_int32 :: s, n
  type(sll_fft_plan), pointer :: fft_plan, fftpack_plan, fftw_plan
  mark => new_time_mark() 
  
  do s=18,20
    n=2**s
    print *,'n=',n
    allocate(copy_data_real(0:n-1))
    allocate(data_real_w(0:n-1))
    allocate(data_real_pack(0:n-1))
    allocate(data_real_sll(0:n-1))

    do i=0,n-1
      phase             = 2.0*sll_pi*real(i)/real(n)
      val               = test_func(phase)
      data_real_w(i)    = val
      data_real_pack(i) = val
      data_real_sll(i)  = val
      copy_data_real(i) = val
    enddo
    
    fft_plan => sll_new_fft(n,FFT_REAL,FFT_NORMALIZE_INVERSE)
    mark => reset_time_mark(mark)
    fft_plan => sll_apply_fft(fft_plan,data_real_sll,FFT_FORWARD)
    time = time_elapsed_since(mark)
    print *, 'SLL_FFT : fft time : ',time
 
    fft_plan => sll_apply_fft(fft_plan,data_real_sll,FFT_INVERSE)
    
    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_sll(i))
    enddo
    print * ,'Averager error: ', acc/n
    print * ,'Max error: ', MAXVAL(ABS(data_real_sll - copy_data_real))


#ifndef _NOFFTW

    fftw_plan => sll_new_fft(n,FFT_REAL,FFTW_MOD + FFT_NORMALIZE_INVERSE)
    mark => reset_time_mark(mark)
    fftw_plan => sll_apply_fft(fftw_plan,data_real_w,FFT_FORWARD)
    time = time_elapsed_since(mark)
    print *, 'FFTW : fft time : ',time

    fftw_plan => sll_apply_fft(fftw_plan,data_real_w,FFT_INVERSE)

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_w(i))
    enddo
    print * ,'Averager error: ', acc/n
    print * ,'Max error: ', MAXVAL(ABS(data_real_w - copy_data_real))


#endif

    fftpack_plan => sll_new_fft(n,FFT_REAL,FFTPACK_MOD + FFT_NORMALIZE_INVERSE)
    mark => reset_time_mark(mark)
    fftpack_plan => sll_apply_fft(fftpack_plan,data_real_pack,FFT_FORWARD)
    time = time_elapsed_since(mark)
    print *, 'FFTPACK : fft time : ',time

    fftpack_plan => sll_apply_fft(fftpack_plan,data_real_pack,FFT_INVERSE)

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_pack(i))
    enddo
    print * ,'Averager error: ', acc/n
    print * ,'Max error: ', MAXVAL(ABS(data_real_pack - copy_data_real))


    deallocate(copy_data_real)
    deallocate(data_real_w)
    deallocate(data_real_pack)
    deallocate(data_real_sll)

    fft_plan => sll_delete_fft(fft_plan)
#ifndef _NOFFTW  
    fftw_plan => sll_delete_fft(fftw_plan)
#endif  
    fftpack_plan => sll_delete_fft(fftpack_plan)  
  enddo
 


contains

  function f(x) result(y)
    sll_int32 :: x
    sll_real64 :: y

    y = real(x,kind=f64)
  end function

  function test_func(x)
    sll_real64, intent(in) :: x
    sll_real64 :: test_func
    test_func = 1.0 + 1.0*cos(x) + 2.0*cos(2.0*x) + 3.0*cos(3.0*x) + &
         4.0*cos(4.0*x) + 5.0*cos(5.0*x) + 6.0*cos(6.0*x) + 7.0*cos(7.0*x) + &
         8.0*cos(8.0*x) + &
         1.0*sin(x) + 2.0*sin(2.0*x) + 3.0*sin(3.0*x) + &
         4.0*sin(4.0*x) + 5.0*sin(5.0*x) + 6.0*sin(6.0*x) + 7.0*sin(7.0*x) + &
         8.0*sin(8.0*x) 
  end function test_func
end program unit_test
