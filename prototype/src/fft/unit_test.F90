program unit_test
  use sll_fft
  use numeric_constants
  use sll_timer

#include "sll_working_precision.h"
  implicit none
  
  sll_int32 :: i, k
  sll_real64, allocatable, dimension(:)  :: copy_data_real, data_real_sll, data_real_w, data_real_pack
  sll_comp64, allocatable, dimension(:)  :: copy_data_comp, data_comp_sll, data_comp_w, data_comp_pack
  type(time_mark), pointer :: mark 
  sll_real64 :: val, phase, acc, acc2
  type(sll_fft_plan), pointer :: fft_plan, fftpack_plan, fftw_plan
  sll_comp64 :: diff,  mode
  sll_int32 :: s, n, t, t_max, s_max, s_min
  sll_real64, allocatable, dimension(:) :: times, times2, times3, accuracy

  mark => new_time_mark() 
 
#define _PRINT
!#define _WRITE
#define _TESTREALONLY
#define _WRITEACCESS

  s_min = 15
  s_max = 20
  t_max = 1

  allocate(times(t_max))
  allocate(times2(t_max))
  allocate(times3(t_max))
  allocate(accuracy(t_max))

#ifdef _PRINT
  print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
#endif

#ifndef _TESTREALONLY
  ! TEST FOR FFTPACK
  open(1,file='fftpackc.txt')
  do s=s_min,s_max
   do t=1,t_max
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_comp(0:n-1))
    allocate(data_comp_pack(0:n-1))

    do i=0,n-1
      phase             = 2.0_f64*sll_pi*real(i,kind=f64)/real(n,kind=f64)
      valc              = test_func_comp(phase)
      data_comp_pack(i) = valc
      copy_data_comp(i) = valc
    enddo

    fftpack_plan => fft_plan_c2c_1d(FFTPACK_MOD,n,data_comp_pack,data_comp_pack,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_c2c_1d(fftpack_plan,data_comp_pack,data_comp_pack)
    times(t) = time_elapsed_since(mark)
    fftpack_plan => fft_plan_c2c_1d(FFTPACK_MOD,n,data_comp_pack,data_comp_pack,FFT_INVERSE)
    call fft_apply_c2c_1d(fftpack_plan,data_comp_pack,data_comp_pack)

    data_comp_pack = data_comp_pack/real(n,kind=f64)

    diff = cmplx(0.0_f64,0.0_f64,kind=f64)
    do i=0,n-1
      diff = diff + abs(copy_data_comp(i) - data_comp_pack(i))
    enddo
    accuracy(t) =  (real(diff) + imag(diff))/(2.0_f64*real(n,kind=f64))

    ! Time access mode
    mark => reset_time_mark(mark)    
    do i=0,n-1
      !k = i!sll_get_index_complex_array(fftpack_plan,i)
      !mode = sll_get_mode_complex_array(fftpack_plan,data_comp_pack,0)
      mode = data_comp_pack(i)
    enddo
    times2(t) = time_elapsed_since(mark)
    print*,times2(t)

    mark => reset_time_mark(mark)    
    do i=0,n-1
      mode = data_comp_pack(i)
    enddo
    times3(t) = time_elapsed_since(mark)
    print*,times3(t)

    deallocate(copy_data_comp)
    deallocate(data_comp_pack)
    call delete(fftpack_plan)
   enddo
#ifdef _PRINT
  print *, 'FFTPACK',SUM(times)/real(t_max,kind=f64),&
                    (SUM(times3)-SUM(times2))/real(t_max,kind=f64),&
                     SUM(accuracy)/real(t_max,kind=f64)
#endif
#ifdef _WRITE
  write(1,*) s,SUM(times)/real(t_max,kind=f64)
#endif
  enddo
  close(1)
  ! END TEST FFTACK



  ! TEST FOR SELALIB FFT
  open(2,file='sllc.txt')
  do s=s_min,s_max
   do t=1,t_max
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_comp(0:n-1))
    allocate(data_comp_sll(0:n-1))

    do i=0,n-1
      phase             = 2.0_f64*sll_pi*real(i,kind=f64)/real(n,kind=f64)
      valc              = test_func_comp(phase)
      data_comp_sll(i)  = valc
      copy_data_comp(i) = valc
    enddo

    fft_plan => fft_plan_c2c_1d(SLLFFT_MOD,n,data_comp_sll,data_comp_sll,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_c2c_1d(fft_plan,data_comp_sll,data_comp_sll)
    times(t) = time_elapsed_since(mark)

    !call CFFTF(copy_data_comp)
    !print *,'==>',copy_data_comp - data_comp_sll

    fft_plan => fft_plan_c2c_1d(SLLFFT_MOD,n,data_comp_sll,data_comp_sll,FFT_INVERSE)
    call fft_apply_c2c_1d(fft_plan,data_comp_sll,data_comp_sll)

    data_comp_sll = data_comp_sll/real(n,kind=f64)

    diff = cmplx(0.0_f64,0.0_f64,kind=f64)
    do i=0,n-1
      diff = diff + abs(copy_data_comp(i) - data_comp_sll(i))
    enddo
    accuracy(t) =  (real(diff) + imag(diff))/(2.0_f64*real(n,kind=f64))
    deallocate(copy_data_comp)
    deallocate(data_comp_sll)
    call delete(fft_plan)
   enddo
#ifdef _PRINT
  print *, ' SLLFFT',SUM(times)/real(t_max,kind=f64),0.d0, SUM(accuracy)/real(t_max,kind=f64)
#endif
#ifdef _WRITE
  write(2,*) s,SUM(times)/real(t_max,kind=f64)
#endif
  enddo
  close(2)
! END SELALIB FFT


#ifndef _NOFFTW
  open(3,file='fftwc.txt')
  do s=s_min,s_max
   do t=1,t_max
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_comp(0:n-1))
    allocate(data_comp_w(0:n-1))

    do i=0,n-1
      phase             = 2.0_f64*sll_pi*real(i,kind=f64)/real(n,kind=f64)
      valc              = test_func_comp(phase)
      data_comp_w(i)    = valc
      copy_data_comp(i) = valc
    enddo

    fftw_plan => fft_plan_c2c_1d(FFTW_MOD,n,data_comp_w,data_comp_w,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_c2c_1d(fftw_plan,data_comp_w,data_comp_w)
    times(t) = time_elapsed_since(mark)

    fftw_plan => fft_plan_c2c_1d(FFTW_MOD,n,data_comp_w,data_comp_w,FFT_INVERSE)
    call fft_apply_c2c_1d(fftw_plan,data_comp_w,data_comp_w)

    data_comp_w = data_comp_w/real(n,kind=f64)

    diff = cmplx(0.0_f64,0.0_f64,kind=f64)
    do i=0,n-1
      diff = diff + abs(copy_data_comp(i) - data_comp_w(i))
    enddo
    accuracy(t) =  (real(diff) + imag(diff))/(2.0_f64*real(n,kind=f64))
    deallocate(copy_data_comp)
    deallocate(data_comp_w)
    call delete(fftw_plan)
   enddo
#ifdef _PRINT
  print *, '   FFTW',SUM(times)/real(t_max,kind=f64),0.d0, SUM(accuracy)/real(t_max,kind=f64)
#endif
#ifdef _WRITE
   write(3,*) s,SUM(times)/real(t_max,kind=f64)
#endif
  enddo
  close(3)
#endif
#endif



#define _NOFFTPACK
#define _NOFFTW







          !!!!!!!!!!!!!!!!!
          ! TEST FOR REAL !
          !!!!!!!!!!!!!!!!!








#ifndef _NOFFTPACK
  ! TEST FOR FFTPACK
  open(4,file='fftpack.txt')
  do s=s_min,s_max
   do t=1,t_max
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_real(0:n-1))
    allocate(data_real_pack(0:n-1))

    do i=0,n-1
      phase             = 2.0_f64*sll_pi*real(i,kind=f64)/real(n,kind=f64)
      val               = test_func(phase)
      data_real_pack(i) = val
      copy_data_real(i) = val
    enddo
    fftpack_plan => fft_plan_r2r_1d(FFTPACK_MOD,n,data_real_pack,data_real_pack,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_r2r_1d(fftpack_plan,data_real_pack,data_real_pack)
    times(t) = time_elapsed_since(mark)
    fftpack_plan => fft_plan_r2r_1d(FFTPACK_MOD,n,data_real_pack,data_real_pack,FFT_INVERSE)
    call fft_apply_r2r_1d(fftpack_plan,data_real_pack,data_real_pack)

    data_real_pack = data_real_pack/real(n,kind=f64)

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_pack(i))
    enddo
    accuracy(t) = acc/real(n,kind=f64)

    ! Time access mode
    mark => reset_time_mark(mark)    
    do i=0,n-1
      !k = sll_get_index_real_array(fftpack_plan,i)
      mode = sll_get_mode_real_array(fftpack_plan,data_real_pack,&
                                     sll_get_index_real_array(fftpack_plan,i))
    enddo
    times2(t) = time_elapsed_since(mark)

    mark => reset_time_mark(mark)
    acc = 0.0_f64    
    !mode k=0
    mode = cmplx(data_real_pack(0),0.0_f64,kind=f64)
    val =  ABS(mode)
    acc = acc + val
    !mode k=1 to k=n-2
    do i=1,n-2
      if(i .eq. (i+1)/2) then !i even
        mode = cmplx(data_real_pack(i-1),-data_real_pack(i),kind=f64)
        val =  ABS(mode)
        acc = acc + val
      else
        mode = cmplx(data_real_pack(i),data_real_pack(i+1),kind=f64)
        val =  ABS(mode)
        acc = acc + val
      endif
    enddo
    !mode k=n/2
    mode = cmplx(data_real_pack(n-1),0.0_f64,kind=f64)
    val =  ABS(mode)
    acc = acc + val
    times3(t) = time_elapsed_since(mark)

    deallocate(copy_data_real)
    deallocate(data_real_pack)
    call delete(fftpack_plan)
   enddo
#ifdef _PRINT
  print *, 'FFTPACK',SUM(times)/real(t_max,kind=f64),&
                    (SUM(times3)-SUM(times2))/real(t_max,kind=f64),&
                     SUM(accuracy)/real(t_max,kind=f64)
#endif
#ifdef _WRITE
  write(4,*) s,SUM(times)/real(t_max,kind=f64)
#endif
  enddo
  close(4)
  ! END TEST FFTACK
#endif










#ifndef _NOFFTSLL
#ifdef _WRITEACCESS
    open(15,file='access_getmode.txt')
    open(16,file='access_direct.txt')
#endif
  ! TEST FOR SELALIB FFT
  open(5,file='sll.txt')
  do s=s_min,s_max
   do t=1,t_max
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_real(0:n-1))
    allocate(data_real_sll(0:n-1))
    
    !Initialization
    do i=0,n-1
      phase             = 2.0_f64*sll_pi*real(i,kind=f64)/real(n,kind=f64)
      val               = test_func(phase)
      data_real_sll(i)  = val
      copy_data_real(i) = val
    enddo

    !FFT FORWARD
    fft_plan => fft_plan_r2r_1d(SLLFFT_MOD,n,data_real_sll,data_real_sll,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_r2r_1d(fft_plan,data_real_sll,data_real_sll)
    times(t) = time_elapsed_since(mark)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Time access mode
    mark => reset_time_mark(mark)
    acc = 0.0_f64    
    do i=0,n-1
      k = sll_get_index_real_array(fft_plan,i)
      mode = sll_get_mode_real_array(fft_plan,data_real_sll,k)
      val =  ABS(mode)
      acc = acc + val
    enddo
    times2(t) = time_elapsed_since(mark)
    !print *,'time with get_mode ==> ', times2
    acc2 = acc

    mark => reset_time_mark(mark)
    acc = 0.0_f64    
    !mode k=0
    mode = cmplx(data_real_sll(0),0.0_f64,kind=f64)
    val =  ABS(mode)
    acc = acc + val
    !mode k=n/2
    mode = cmplx(data_real_sll(1),0.0_f64,kind=f64)
    val =  ABS(mode)
    acc = acc + val
    !mode k=1 to k= n-1
    do i=1,n/2-1
        mode = cmplx(data_real_sll(2*i),data_real_sll(2*i+1),kind=f64)
        val =  ABS(mode)
        acc = acc + val
        mode = cmplx(data_real_sll(2*i),-data_real_sll(2*i+1),kind=f64)
        val =  ABS(mode)
        acc = acc + val
    enddo
    times3(t) = time_elapsed_since(mark)
    !print *,'time with direct access ==> ', times3

    if( acc .ne. acc2 ) &
      print*,'prob=',abs(acc-acc2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !FFT INVERSE
    call delete(fft_plan)
    fft_plan => fft_plan_r2r_1d(SLLFFT_MOD,n,data_real_sll,data_real_sll,FFT_INVERSE,FFT_NORMALIZE_INVERSE)
    call fft_apply_r2r_1d(fft_plan,data_real_sll,data_real_sll)
    !Normalization
    !data_real_sll = 2.0_f64*data_real_sll/real(n,kind=f64)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !ACCURATE
    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_sll(i))
    enddo
    accuracy(t) =  acc/real(n,kind=f64)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    deallocate(copy_data_real)
    deallocate(data_real_sll)
    call delete(fft_plan)
   enddo
#ifdef _PRINT
  print *, ' SLLFFT',SUM(times)/real(t_max,kind=f64),&
                    (SUM(times3)-SUM(times2))/real(t_max,kind=f64),&
                     SUM(accuracy)/real(t_max,kind=f64)
#endif
#ifdef _WRITEACCESS
    write(15,*) n, SUM(times2)/real(t_max,kind=f64)
    write(16,*) n, SUM(times3)/real(t_max,kind=f64)
#endif
#ifdef _WRITE
  write(5,*) s,SUM(times)/real(t_max,kind=f64)
#endif
  enddo
  close(5)
! END SELALIB FFT
#endif
#ifdef _WRITEACCESS
    close(15)
    close(16)
#endif



#ifndef _NOFFTW
  open(7,file='fftw.txt')
  do s=s_min,s_max
   do t=1,t_max
    n=2**s
    !print *,'n=',n
    !print *, 'LIB      TIME EXECUTION            TIME ACCESS MODE          AVERAGE ERROR'
    allocate(copy_data_real(0:n-1))
    allocate(data_real_w(0:n-1))
    do i=0,n-1
      phase             = 2.0_f64*sll_pi*real(i,kind=f64)/real(n,kind=f64)
      val               = test_func(phase)
      data_real_w(i)    = val
      copy_data_real(i) = val
    enddo

    fftw_plan => fft_plan_r2r_1d(FFTW_MOD,n,data_real_w,data_real_w,FFT_FORWARD)
    mark => reset_time_mark(mark)
    call fft_apply_r2r_1d(fftw_plan,data_real_w,data_real_w)
    times(t) = time_elapsed_since(mark)

    fftw_plan => fft_plan_r2r_1d(FFTW_MOD,n,data_real_w,data_real_w,FFT_INVERSE)
    call fft_apply_r2r_1d(fftw_plan,data_real_w,data_real_w)

    data_real_w = data_real_w/real(n,kind=f64)

    acc = 0.0_f64
    do i=0,n-1
      acc = acc + abs(copy_data_real(i) - data_real_w(i))
    enddo
    accuracy(t) =  acc/real(n,kind=f64)

    ! Time access mode
    mark => reset_time_mark(mark)    
    do i=0,n-1
      !k = sll_get_index_real_array(fftw_plan,i)
      mode = sll_get_mode_real_array(fftw_plan,data_real_w,&
                                     sll_get_index_real_array(fftw_plan,i))
    enddo
    times2(t) = time_elapsed_since(mark)

    mark => reset_time_mark(mark)
    acc = 0.0_f64    
    !mode k=0
    mode = cmplx(data_real_pack(0),0.0_f64,kind=f64)
    val =  ABS(mode)
    acc = acc + val
    !mode k=1 to k=n/2-1
    do i=1,n/2-1
      mode = cmplx(data_real_w(i) , data_real_w(n-i),kind=f64)
      val =  ABS(mode)
      acc = acc + val
    enddo
    !mode k=n/2
    mode = cmplx(data_real_pack(n-1),0.0_f64,kind=f64)
    val =  ABS(mode)
    acc = acc + val
    !mode k=n/2+1 to k=n-1
    do i=n/2+1,n-1
      mode = cmplx(data_real_w(n-i) , -data_real_w(i),kind=f64)
      val =  ABS(mode)
      acc = acc + val
    enddo
    times3(t) = time_elapsed_since(mark)

    deallocate(copy_data_real)
    deallocate(data_real_w)
    call delete(fftw_plan)
   enddo
#ifdef _PRINT
  print *, '   FFTW',SUM(times)/real(t_max,kind=f64),&
                    (SUM(times3)-SUM(times2))/real(t_max,kind=f64),&
                     SUM(accuracy)/real(t_max,kind=f64)
#endif
#ifdef _WRITE
  write(7,*) s,SUM(times)/real(t_max,kind=f64)
#endif
  enddo
  close(7)
#endif

  deallocate(accuracy)
  deallocate(times)
  deallocate(times2)
  deallocate(times3)

#ifdef _PRINT
#undef _PRINT
#endif
#ifdef _WRITE
#undef _WRITE
#endif

contains

  function f(x) result(y)
    sll_int32 :: x
    sll_comp64 :: y

    y = cmplx(real(x,kind=f64),0.0_f64,kind=f64)
  end function

  function g(x) result(y)
    sll_real64 :: x
    sll_real64 :: y

    y = x
  end function

  function test_func(x)
    sll_real64, intent(in) :: x
    sll_real64 :: test_func
    test_func = 1.0_f64 + 1.0_f64*cos(x) + 2.0_f64*cos(2.0_f64*x) + 3.0_f64*cos(3.0_f64*x) + &
         4.0_f64*cos(4.0_f64*x) + 5.0_f64*cos(5.0_f64*x) + 6.0_f64*cos(6.0_f64*x) + 7.0_f64*cos(7.0_f64*x) + &
         8.0_f64*cos(8.0_f64*x) + &
         1.0_f64*sin(x) + 2.0_f64*sin(2.0_f64*x) + 3.0_f64*sin(3.0_f64*x) + &
         4.0_f64*sin(4.0_f64*x) + 5.0_f64*sin(5.0_f64*x) + 6.0_f64*sin(6.0_f64*x) + 7.0_f64*sin(7.0_f64*x) + &
         8.0_f64*sin(8.0_f64*x) 
  end function test_func

  function test_func_comp(x)
    sll_real64, intent(in) :: x
    sll_comp64 :: test_func_comp
    test_func_comp = cmplx(1.0_f64 + 1.0_f64*cos(x) + 2.0_f64*cos(2.0_f64*x) + 3.0_f64*cos(3.0_f64*x) + &
         4.0_f64*cos(4.0_f64*x) + 5.0_f64*cos(5.0_f64*x) + 6.0_f64*cos(6.0_f64*x) + 7.0_f64*cos(7.0_f64*x) + &
         8.0_f64*cos(8.0_f64*x) ,&
         1.0_f64*sin(x) + 2.0_f64*sin(2.0_f64*x) + 3.0_f64*sin(3.0_f64*x) + &
         4.0_f64*sin(4.0_f64*x) + 5.0_f64*sin(5.0_f64*x) + 6.0_f64*sin(6.0_f64*x) + 7.0_f64*sin(7.0_f64*x) + &
         8.0_f64*sin(8.0_f64*x) , kind=f64)
  end function test_func_comp

  function getmode(plan,data,k) result(mode)
    sll_comp64                                :: mode
    sll_int32, intent(in)                     :: k
    type(sll_fft_plan) , pointer              :: plan
    sll_real64, dimension(0:) , intent(in)    :: data
    sll_int32                                 :: n_2, n
   
    n = plan%N
    n_2 = ishft(n,-1)

      if( k .eq. 0 ) then
        mode = complex(data(0),0.0_f64)
      else if( k .eq. n_2 ) then
        mode = complex(data(1),0.0_f64)
      else if( k .gt. n_2 ) then
        mode = complex( data(2*(n-k)) , -data(2*(n-k)+1) )
      else
        mode = complex( data(2*k) , data(2*k+1) )
      endif
      return
  end function

end program unit_test
