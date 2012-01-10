program unit_test
  use sll_fft
  use numeric_constants
  use sll_timer

#include "sll_working_precision.h"
  implicit none
  
  
  sll_int32 :: i
  sll_real64, allocatable, dimension(:)  :: dat
  sll_real64, allocatable, dimension(:)  :: c
  sll_real64, allocatable, dimension(:)  :: in
  sll_comp64, allocatable, dimension(:)  :: diff
  type(time_mark), pointer :: mark 
  sll_real64 :: time
  sll_int32 :: s, n
  type(sll_fft_plan), pointer :: fft, fft3
  type(sll_fft_plan), pointer :: fft2  
  mark => new_time_mark() 
  
  do s=18,22
    n=2**s
    print *,'n=',n
    allocate(dat(0:n-1))
    allocate(in(0:n-1))
    allocate(c(0:n-1))
    allocate(diff(0:n-1))

    do i=0,n-1
      dat(i) = f(i)
      in(i) = f(i)
      c(i) = f(i)
    enddo

    fft => sll_new_fft(n,FFT_REAL)
    mark => reset_time_mark(mark)
    fft => sll_apply_fft(fft,dat,FFT_FORWARD)
    time = time_elapsed_since(mark)
    print *, 'SLL_FFT : fft time : ',time

#ifndef _NOFFTW
    fft2 => sll_new_fft(n,FFT_REAL,FFTW_MOD)
    mark => reset_time_mark(mark)
    fft2 => sll_apply_fft(fft2,in,FFT_FORWARD)
    time = time_elapsed_since(mark)
    print *, 'FFTW : fft time : ',time
#endif

    fft3 => sll_new_fft(n,FFT_REAL,FFTPACK_MOD)
    mark => reset_time_mark(mark)
    fft3 => sll_apply_fft(fft3,c,FFT_FORWARD)
    time = time_elapsed_since(mark)
    print *, 'FFTPACK : fft time : ',time

    !open(4,file='dat.txt')
    !do i=0,n-1
    !  write(4,*) dat(i) , c(i)
    !enddo  
    !close(4)

    do i=0,n-1
      diff(i) = abs(sll_get_mode(i,fft,dat) - sll_get_mode(i,fft2,in))
    enddo
    print *, MAXVAL( real(diff) ), MAXVAL( imag(diff) )

    deallocate(dat)
    deallocate(diff)
    deallocate(in)
    deallocate(c)
    fft => sll_delete_fft(fft)
#ifndef _NOFFTW  
    fft2 => sll_delete_fft(fft2)
#endif  
    fft3 => sll_delete_fft(fft3)  
  enddo
 


contains
! for a real input we have a complex output of the form
!      (   r_0 ,   0    )         # 0
!      (   r_1 ,   i_1  )         # 1
!      (   r_2 ,   i_2  )         # 2
!         .         .
!         .         .
!         .         .
! (  r_{n/2-1} ,  i_{n/2-1}  )
! (  r_n/2     ,      0      )    # n/2
! (  r_{n/2-1} ,  -i_{n/2-1} )
! (  r_{n/2-2} ,  -i_{n/2-2} )
!        .           .
!        .           .
!        .           .
!     (   r_2  ,   -i_2 )
!     (   r_1  ,   -i_1 )         # n/2 + n/2 - 1 = n - 1
!
!
! the real output in FFTW have the form
! [r_0,r_1,....,r_n/2,i_(n/2 - 1),...,i_1]
!
! the real output in SLL_FFT have the form
! [r_0,r_n/2,r_1,i_1,r_2,i_2,..,i_(n/2 - 1),i_(n/2 - 1)]
!
! the real output in FFTPACK have the form
! [r_0,,r_1,i_1,r_2,i_2,..,i_(n/2 - 1),i_(n/2 - 1),r_n/2]

  function switch_halfcomplex_sll_to_fftw(n,t)
    sll_real64, dimension(0:n-1)             :: switch_halfcomplex_sll_to_fftw
    sll_int32, intent(in)                    :: n
    sll_real64, dimension(0:n-1), intent(in) :: t
    sll_int32 :: i,j

    switch_halfcomplex_sll_to_fftw(0) = t(0)
    switch_halfcomplex_sll_to_fftw(n/2) = t(1) 
        
    j=0
    do i=2,n-2,2
      j=j+1
      switch_halfcomplex_sll_to_fftw(j) = t(i)
      switch_halfcomplex_sll_to_fftw(n-j) = t(i+1)
    enddo
  end function switch_halfcomplex_sll_to_fftw

  function switch_halfcomplex_fftpack_to_fftw(n,t)
    sll_real64, dimension(0:n-1)             :: switch_halfcomplex_fftpack_to_fftw
    sll_int32, intent(in)                    :: n
    sll_real64, dimension(0:n-1), intent(in) :: t
    sll_int32 :: i,j

    switch_halfcomplex_fftpack_to_fftw(0) = t(0)
    switch_halfcomplex_fftpack_to_fftw(n/2) = t(n-1) 
        
    j=0
    do i=1,n-2,2
      j=j+1
      switch_halfcomplex_fftpack_to_fftw(j) = t(i)
      switch_halfcomplex_fftpack_to_fftw(n-j) = t(i+1)
    enddo
  end function switch_halfcomplex_fftpack_to_fftw

  function f(x) result(y)
    sll_int32 :: x
    sll_real64 :: y

    y = cos(real(x,kind=f64))
  end function
end program unit_test
