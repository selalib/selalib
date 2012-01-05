program unit_test
  use sll_fft
  use numeric_constants
  use sll_timer

#include "sll_working_precision.h"
  implicit none
  
  
  sll_int32 :: i
  sll_real64, allocatable, dimension(:)  :: dat, c
  !sll_real64, allocatable, dimension(:)  :: in
  type(time_mark), pointer :: mark 
  sll_real64 :: time
  sll_int32 :: s, n
  type(sll_fft_plan), pointer :: fft, fft3
 ! type(sll_fft_plan), pointer :: fft2  
  mark => new_time_mark() 

  do s=20,25
    n=2**s
    print *,'n=',n
    allocate(dat(0:n-1))
    !allocate(in(0:n-1))
    allocate(c(0:n-1))

    dat(0) = 1.d0
    c(0) = 1.d0
    !in(0) = 1.d0
    do i=1,n-1
      dat(i) = 0.d0 !real(i,kind=f64)
      c(i) = dat(i)
      !in(i) = dat(i)
    enddo

    mark => reset_time_mark(mark)
    fft => sll_new_fft(n,FFT_REAL)
    time = time_elapsed_since(mark)
    print *, 'SLL_FFT : initializing time : ',time
    mark => reset_time_mark(mark)
    fft => sll_apply_fft(fft,dat,FFT_FORWARD)
    time = time_elapsed_since(mark)
    print *, 'SLL_FFT : fft time : ',time

    !mark => start_time_mark(mark)
    !fft2 => sll_new_fft(n,FFT_REAL,FFTW_MOD)
    !time = time_elapsed_since(mark)
    !print *, 'FFTW : initializing time : ',time
    !mark => reset_time_mark(mark)
    !fft2 => sll_apply_fft(fft2,in,FFT_FORWARD)
    !time = time_elapsed_since(mark)
    !print *, 'FFTW : fft time : ',time

    mark => reset_time_mark(mark)
    fft3 => sll_new_fft(n,FFT_REAL,FFTPACK_MOD)
    time = time_elapsed_since(mark)
    print *, 'FFTPACK : initializing time : ',time
    mark => reset_time_mark(mark)
    fft3 => sll_apply_fft(fft3,c,FFT_FORWARD)
    time = time_elapsed_since(mark)
    print *, 'FFTPACK : fft time : ',time

    dat = switch_halfcomplex_sll_to_fftw(n,dat)
    c = switch_halfcomplex_fftpack_to_fftw(n,c)

    !open(4,file='dat.txt')
    !do i=0,n-1
    !  write(4,*) dat(i) , in(i), c(i)!, r(i)
    !enddo  
    !close(4)

    !print *, MAXVAL(ABS(2*in-dat-c))
    print *, MAXVAL(ABS(dat-c))
    deallocate(dat)
    !deallocate(in)
    deallocate(c)
    fft => sll_delete_fft(fft)  
    !fft2 => sll_delete_fft(fft2)  
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
end program unit_test
