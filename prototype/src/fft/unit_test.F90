program unit_test
  use sll_fft
  use numeric_constants
  use sll_timer

#include "sll_working_precision.h"
  implicit none
  
  sll_int32 :: i
  sll_real64, allocatable, dimension(:)  :: copy_data_real, data_real_sll, data_real_w, data_real_pack
  type(time_mark), pointer :: mark 
  sll_real64 :: val, phase, acc 
  sll_int32 :: s, n, t
  sll_real64, dimension(10) :: times
  type(sll_fft_plan), pointer :: fft_plan, fftpack_plan, fftw_plan
  mark => new_time_mark() 
    

  ! TEST FOR FFTPACK
  open(1,file='fftpack.txt')
  do s=10,20
   do t=1,10
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_real(0:n-1))
    allocate(data_real_pack(0:n-1))

    do i=0,n-1
      phase             = 2.0*sll_pi*real(i)/real(n)
      val               = g(phase)
      data_real_pack(i) = val
      copy_data_real(i) = val
    enddo

    fftpack_plan => fft_plan_r2r_1d(FFTPACK_MOD,n,data_real_pack,data_real_pack,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_r2r_1d(fftpack_plan,data_real_pack,data_real_pack)
    times(t) = time_elapsed_since(mark)
    fftpack_plan => fft_plan_r2r_1d(FFTPACK_MOD,n,data_real_pack,data_real_pack,FFT_INVERSE)
    call fft_apply_r2r_1d(fftpack_plan,data_real_pack,data_real_pack)

    data_real_pack = data_real_pack/n

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_pack(i))
    enddo
    !print *, 'FFTPACK',time,time, acc/n
 
    deallocate(copy_data_real)
    deallocate(data_real_pack)
    call fft_delete_plan(fftpack_plan)
   enddo
  write(1,*) s,SUM(times)/10.0d0
  enddo
  close(1)
  ! END TEST FFTACK



  ! TEST FOR SELALIB FFT
  open(2,file='sll.txt')
  do s=10,20
   do t=1,10
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_real(0:n-1))
    allocate(data_real_sll(0:n-1))

    do i=0,n-1
      phase             = 2.0*sll_pi*real(i)/real(n)
      val               = g(phase)
      data_real_sll(i)  = val
      copy_data_real(i) = val
    enddo

    fft_plan => fft_plan_r2r_1d(SLLFFT_MOD,n,data_real_sll,data_real_sll,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_r2r_1d(fft_plan,data_real_sll,data_real_sll)
    times(t) = time_elapsed_since(mark)

    fft_plan => fft_plan_r2r_1d(SLLFFT_MOD,n,data_real_sll,data_real_sll,FFT_INVERSE)
    call fft_apply_r2r_1d(fft_plan,data_real_sll,data_real_sll)

    data_real_sll = 2.0_f64*data_real_sll/n

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_sll(i))
    enddo
    !print *, ' SLLFFT',time,time, acc/n

    deallocate(copy_data_real)
    deallocate(data_real_sll)
    call fft_delete_plan(fft_plan)
   enddo
   write(2,*) s,SUM(times)/10.0d0
  enddo
  close(2)
! END SELALIB FFT

#ifndef _NOFFTW
  open(3,file='fftw.txt')
  do s=10,20
   do t=1,10
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_real(0:n-1))
    allocate(data_real_w(0:n-1))

    do i=0,n-1
      phase             = 2.0*sll_pi*real(i)/real(n)
      val               = g(phase)
      data_real_w(i)    = val
      copy_data_real(i) = val
    enddo

    fftw_plan => fft_plan_r2r_1d(FFTW_MOD,n,data_real_w,data_real_w,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_r2r_1d(fftw_plan,data_real_w,data_real_w)
    times(t) = time_elapsed_since(mark)

    fftw_plan => fft_plan_r2r_1d(FFTW_MOD,n,data_real_w,data_real_w,FFT_INVERSE)
    call fft_apply_r2r_1d(fftw_plan,data_real_w,data_real_w)

    data_real_w = data_real_w/n

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_w(i))
    enddo
    !print *, '   FFTW',time,time, acc/n

    deallocate(copy_data_real)
    deallocate(data_real_w)
    call fft_delete_plan(fftw_plan)
   enddo
   write(3,*) s,SUM(times)/10.0d0
  enddo
  close(3)
#endif


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
