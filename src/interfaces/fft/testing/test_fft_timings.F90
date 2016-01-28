! Test program to test performance of FFTW interface type-bound versus non-type bound.
! Test for FFTW-versions on Intel Ivy Bridge notebook processor 3.0 GHz with gfortran 4.8.4
! Computing times (for initialize, exexute, delete): min / average / max (out of 10 runs)
! Non-type bound version: 0.3082 / 0.3136 / 0.3178
! Type bound version: 0.3132 / 0.3165 / 0.3233
! Type bound version but called non-type bound: 0.3132 / 0.3164 / 0.3218
! Performance decrease for type bound: 1.6% / 0.9% /  1.7%
! AUTHOR: Katharina Kormann, IPP
! Date: January 25, 2016

! Timings for just apply function:
! Non-type bound version: 0.0088 / 0.0119 / 0.0170
! Type bound version but called non-type bound: 0.0109 / 0.0135 / 0.0160
! Performance decrease for type bound: 24% / 13% /-6%

! Timing with 10000000 instead of 100000 calls:
! 2.5% overhead on minimum and 3.3 on average
! Using select type within subroutine brings overhead down to 
! 1.1% overhead on minimum and 1.5 on average
! Same on hydra (with intel 14): 9.9 % overhead on minimum and 10.2 % on average
! Type:  0.8065 / 0.8254 / 0.8456
! Class: 0.8863 / 0.9098 / 0.9353

program test_fft_timing
#include "sll_working_precision.h"

  use sll_m_fft

  use sll_m_timer

  implicit none

  type(sll_t_fft) :: p

  sll_int32, parameter :: n = 16

  sll_comp64 :: data_in(n)
  sll_comp64 :: data_out(n)

  sll_int32 :: rnd_seed_size
  sll_int32, allocatable :: rnd_seed(:) !< Random seed.

  sll_real64 :: time
  type(sll_t_time_mark) :: t0

  sll_int32 :: j

 ! Set some random seed
  call random_seed (size=rnd_seed_size)
  allocate(rnd_seed(rnd_seed_size))
  do j=1, rnd_seed_size
     rnd_seed(j) = 100+15*j
  end do
  call random_seed (put=rnd_seed)

  do j=1,n
     call RANDOM_COMPLEX(data_in(j))
  end do
  

  call sll_s_fft_init_c2c_1d(p, n, data_in, data_out, &
       sll_p_fft_forward)
  call sll_s_fft_exec_c2c_1d(p, data_in, data_out)
  call sll_s_set_time_mark(t0)
  do j=1,10000000
     call sll_s_fft_exec_c2c_1d(p, data_in, data_out)
     
     !print*, 'a'
  end do
  
  time = sll_f_time_elapsed_since(t0)
  call sll_s_fft_free(p)
 ! print*, 'Computing time for fft c2c functions', time
  print*, time
  
!!$  call sll_s_set_time_mark(t0)
!!$  do j=1,100000
!!$     call p%initialize( n, data_in, data_out, &
!!$          sll_p_fft_forward)
!!$     call p%execute( data_in, data_out)
!!$     
!!$     call p%delete()
!!$     !print*, 'a'
!!$  end do
!!$  time = sll_f_time_elapsed_since(t0)
!!$  print*, 'Computing time for fft c2c functions (type bound)', time




contains

  SUBROUTINE init_random_seed()
    INTEGER                            :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
  END SUBROUTINE

  SUBROUTINE RANDOM_COMPLEX(c)
    sll_comp64, intent(out) :: c
    sll_real64 :: realpart, imagpart
   
    CALL init_random_seed() 
    CALL RANDOM_NUMBER(realpart)
    CALL RANDOM_NUMBER(imagpart)
    c = CMPLX(realpart,imagpart,kind=f64)
  END SUBROUTINE

end program test_fft_timing
