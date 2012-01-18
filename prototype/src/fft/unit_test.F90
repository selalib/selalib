program unit_test
  use sll_fft
  use numeric_constants
  use sll_timer

#include "sll_working_precision.h"
  implicit none
 
!#define _NOFFTW 
  
  sll_int32 :: i
  sll_real64, allocatable, dimension(:)  :: copy_data_real, data_real_sll, data_real_w, data_real_pack
  type(time_mark), pointer :: mark 
  sll_real64 :: time, val, phase, acc, time2, time3
  sll_int32 :: s, n
  type(sll_fft_plan), pointer :: fft_plan, fftpack_plan, fftw_plan
  type(C_PTR) :: p, plan
  real(C_DOUBLE), pointer :: a(:)
  sll_real64, dimension(2**22) :: arr
  sll_int32 :: ind
  sll_comp64 :: mode
  
  p = fftw_alloc_real(int(2**20,C_SIZE_T))
  call c_f_pointer(p,a,[2**20])

  mark => new_time_mark() 
   
  do s=20,20
    n=2**s
    print *,'n=',n
    print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_real(0:n-1))
    allocate(data_real_w(0:n-1))
    allocate(data_real_pack(0:n-1))
    allocate(data_real_sll(0:n-1))

    do i=0,n-1
      phase             = 2.0*sll_pi*real(i)/real(n)
      val               = g(phase)
      data_real_w(i)    = val
      data_real_pack(i) = val
      data_real_sll(i)  = val
      a(i+1)            = val
      arr(i+1)          = val
      copy_data_real(i) = val
    enddo
    
    fft_plan => sll_new_fft(n,FFT_REAL,FFT_NORMALIZE_INVERSE)
    mark => reset_time_mark(mark)
    fft_plan => sll_apply_fft(fft_plan,data_real_sll,FFT_FORWARD)
    time = time_elapsed_since(mark)
    !print *, 'SLL_FFT : fft time : ',time
 
    fft_plan => sll_apply_fft(fft_plan,data_real_sll,FFT_INVERSE)

    mark => reset_time_mark(mark)    
    do i=0,n-1
      ind = sll_get_index(i,fft_plan)
      mode = sll_get_mode(i,fft_plan,data_real_sll)
    enddo 
    time2 = time_elapsed_since(mark)
    !print *, 'SLL_FFT : access time : ',time

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_sll(i))
    enddo
    !print * ,'Averager error: ', acc/n
    !print * ,'Max error: ', MAXVAL(ABS(data_real_sll - copy_data_real))
    print *, '    SLL',time,time2,acc/n


#ifndef _NOFFTW
    fftw_plan => sll_new_fft(n,FFT_REAL,FFTW_MOD + FFT_NORMALIZE_INVERSE)
    mark => reset_time_mark(mark)
    fftw_plan => sll_apply_fft(fftw_plan,data_real_w,FFT_FORWARD)
    time = time_elapsed_since(mark)
    !print *, FFTW : fft time : ',time

    fftw_plan => sll_apply_fft(fftw_plan,data_real_w,FFT_INVERSE)

    mark => reset_time_mark(mark)    
    do i=0,n-1
      ind = sll_get_index(i,fftw_plan)
      mode = sll_get_mode(i,fftw_plan,data_real_w)
    enddo 
    time2 = time_elapsed_since(mark)
    !print *, 'FFTW_FFT : access time : ',time


    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_w(i))
    enddo
    !print * ,'Averager error: ', acc/n
    !print * ,'Max error: ', MAXVAL(ABS(data_real_w - copy_data_real))
    print *, '   FFTW',time,time2,acc/n
#endif


#define _NOFFTW
#ifndef _NOFFTW
    plan = fftw_plan_r2r_1d(n,arr,arr,FFTW_R2HC,FFTW_ESTIMATE + FFTW_UNALIGNED)
    mark => reset_time_mark(mark)
    call fftw_execute_r2r( plan , arr , arr )
    time = time_elapsed_since(mark)
    print *, 'FFTW : fft time : ',time

    plan = fftw_plan_r2r_1d(n,arr,arr,FFTW_HC2R,FFTW_ESTIMATE)
    call fftw_execute_r2r( plan , arr , arr )

    acc = 0.0_f64
    do i=0,n-1
      arr(i+1) = arr(i+1)/n
      acc = acc + abs(copy_data_real(i) - arr(i+1))
    enddo
    print * ,'Averager error: ', acc/n
    !print * ,'Max error: ', MAXVAL(ABS(a - copy_data_real))
#endif
#undef _NOFFTW

#ifndef _NOFFTW
    plan = fftw_plan_r2r_1d(n,a,a,FFTW_R2HC,FFTW_ESTIMATE)
    mark => reset_time_mark(mark)
    call fftw_execute_r2r( plan , a , a )
    time = time_elapsed_since(mark)
    !print *, 'FFTW : fft time : ',time

    plan = fftw_plan_r2r_1d(n,a,a,FFTW_HC2R,FFTW_ESTIMATE)
    call fftw_execute_r2r( plan , a , a )

    mark => reset_time_mark(mark)    
    do i=0,n-1
      ind = sll_get_index(i,fftw_plan)
      mode = sll_get_mode(i,fftw_plan,a)
    enddo
    time2 = time_elapsed_since(mark)


    acc = 0.0_f64
    do i=0,n-1
      a(i+1) = a(i+1)/n
      acc = acc + abs(copy_data_real(i) - a(i+1))
    enddo
    !print * ,'Averager error: ', acc/n
    !print * ,'Max error: ', MAXVAL(ABS(a - copy_data_real))
    print *, 'FFTWOPT',time,time2,acc/n
#endif

    fftpack_plan => sll_new_fft(n,FFT_REAL,FFTPACK_MOD + FFT_NORMALIZE_INVERSE)
    mark => reset_time_mark(mark)
    fftpack_plan => sll_apply_fft(fftpack_plan,data_real_pack,FFT_FORWARD)
    time = time_elapsed_since(mark)
    !print *, 'FFTPACK : fft time : ',time

    fftpack_plan => sll_apply_fft(fftpack_plan,data_real_pack,FFT_INVERSE)

    mark => reset_time_mark(mark)    
    do i=0,n-1
      ind = sll_get_index(i,fftpack_plan)
      mode = sll_get_mode(i,fftpack_plan,data_real_pack)
    enddo 
    time2 = time_elapsed_since(mark)
    !print *, 'FFTPACK_FFT : access time : ',time

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_pack(i))
    enddo
    !print * ,'Averager error: ', acc/n
    !print * ,'Max error: ', MAXVAL(ABS(data_real_pack - copy_data_real))
    print *, 'FFTPACK',time,time2,acc/n

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

  function g(x) result(y)
    sll_real64 :: x
    sll_real64 :: y

    y = x
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
