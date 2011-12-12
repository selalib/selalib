program fft_test
  use sll_fft
  use numeric_constants
#include "sll_working_precision.h"
  implicit none

#define CREAL(array, i) array(2*(i))
#define CIMAG(array, i) array(2*(i)+1)

#define NC 4096

#define PRINT_ARRAYS 0

  sll_comp64, dimension(:), allocatable :: dat
  sll_comp64, dimension(:), allocatable :: dat_copy
  sll_comp64, dimension(:), allocatable :: tst
  sll_comp64, dimension(:), allocatable :: data_complex
  sll_comp64, dimension(:), allocatable :: data_complex_copy
  sll_comp64, dimension(:), allocatable :: data_complex_n
  sll_comp64, dimension(:), allocatable :: data_complex_n_copy
  sll_real64, dimension(:), allocatable :: data_real
  sll_real64, dimension(:), allocatable :: data_real_copy
  sll_real64, dimension(:), allocatable :: data_real_copy2
  sll_real64, dimension(:), allocatable :: nr_data
  sll_real64, dimension(:), allocatable :: nr_data_copy
  ! twiddles are named according to the size of the FFT that they are
  ! meant to be used with.
  sll_comp64, dimension(:), allocatable :: twiddles_c_n_2
  sll_comp64, dimension(:), allocatable :: twiddles_c_n
  sll_real64, dimension(:), allocatable :: twiddles_r_n_2
  sll_real64, dimension(:), allocatable :: twiddles_r_n
  sll_real64, dimension(:), allocatable :: tst2
  sll_real64                            :: phase
  sll_comp64                            :: tmp_c
  sll_real64                            :: tmp_r
  sll_real64                            :: tmp_i
  sll_real64                            :: acc
  sll_real64                            :: val
  integer                            :: i

#if TEST_TIME
  sll_real32       :: tarray(1:2)
  sll_real32       :: tresult
#endif
  ! caution: size of twiddles is half the size of the array whose fft is       
  ! needed. This is a source of errors.

  print *, 'Array size for tests: ', NC
  allocate(tst(0:NC/2))
  allocate(tst2(0:NC))
  allocate(dat(0:NC-1))
  allocate(dat_copy(0:NC-1))
  allocate(data_real(0:NC-1))
  allocate(data_real_copy(0:NC-1))
  allocate(data_real_copy2(0:NC-1))
  allocate(nr_data(0:NC-1))
  allocate(nr_data_copy(0:NC-1))
  allocate(data_complex(0:NC/2-1))
  allocate(data_complex_copy(0:NC/2-1))
  allocate(data_complex_n(0:NC-1))
  allocate(data_complex_n_copy(0:NC-1))
  allocate(twiddles_c_n(0:NC/2-1))
  allocate(twiddles_c_n_2(0:NC/4-1))
  allocate(twiddles_r_n_2(0:NC/2-1))! these are real arrays storing complex nums
  allocate(twiddles_r_n(0:NC-1)) ! real array storing complex nums
  call compute_twiddles(NC/2,twiddles_c_n_2(0:NC/2-1))
  call compute_twiddles_real_array( NC/2, twiddles_r_n_2(0:NC-1) )
  call compute_twiddles_real_array( NC, twiddles_r_n(0:NC-1) ) 
  do i=0,NC-1
     phase         = 2.0*sll_pi*real(i)/real(NC)
     dat(i)        = test_func(phase)
     data_real(i)  = test_func(phase)
     nr_data(i)    = test_func(phase)
  end do

  print *, '0--------------------------------------------------------------'
  print *, 'test of bit-reversing subroutines. Note: these only work when the arrays are of size 16' ! ... so don't dynamically allocate them...
  tst(:) = (/(1000,0), (1001,0), (1010,0), (1011,0), (1100,0), (1101,0), (1110,0), (1111,0)/)
#if PRINT_ARRAYS
  print *, tst(:)
#endif
  tst2(:) = (/ 1000, 0, 1001, 0, 1010, 0, 1011, 0, 1100, 0, 1101, 0, 1110, 0, 1111, 0 /)

  call bit_reverse_complex(NC/2,tst)
#if PRINT_ARRAYS
  print *, 'bit reversed:'
  do i=0, size(tst)-1
     print *, tst(i)
  end do
  print *, 'bit reversing in pairs test:'
  print *, 'original:'
  do i=0, size(tst2)/2-1
     write (*,'(f10.0, f10.0)') tst2(2*i), tst2(2*i+1)
  end do
  print *, 'bit-reversed:'
  call bit_reverse_in_pairs(size(tst2)/2,tst2)
  do i=0, size(tst2)/2-1
     write (*,'(f10.0, f10.0)') tst2(2*i), tst2(2*i+1)
  end do
#endif

  print *, '1--------------------------------------------------------------'
  print *, 'Test complex FFT'
  do i=0, NC/2-1
     data_complex(i) = transfer(data_real(2*i:),data_complex(0))
  end do
  data_complex_copy(:) = data_complex(:)

#if PRINT_ARRAYS
  print *, 'complex data: '
  do i=0, NC/2-1
     print *, data_complex(i)
  end do
#endif

  call bit_reverse( NC/4, twiddles_c_n_2(0:(NC/4-1)) )
  call fft_dit_nr_aux( data_complex(0:NC/2-1), NC/2, twiddles_c_n_2(0:NC/4-1), &
                       0, FFT_FORWARD )
#if PRINT_ARRAYS
  print *, 'transformed complex array: '
  do i=0, NC/2-1
     print *, data_complex(i)
  end do
#endif

  print *, 'applying inverse transform: '
  call bit_reverse( NC/4, twiddles_c_n_2(0:(NC/4-1)) )
  call fft_dit_rn_aux( data_complex(0:NC/2-1), NC/2, twiddles_c_n_2(0:NC/4-1), &
       1, FFT_INVERSE )
  data_complex(:)      = data_complex(:)*1.0/real(NC/2)
  data_complex_copy(:) = data_complex_copy(:) - data_complex(:)
  tmp_c = (0.0,0.0)
  do i=0,size(data_complex_copy)-1
     tmp_r = tmp_r + abs( real(data_complex_copy(i)))
     tmp_i = tmp_i + abs(aimag(data_complex_copy(i)))
  end do
  tmp_r = tmp_r/size(data_complex_copy)
  tmp_i = tmp_i/size(data_complex_copy)
  print *, 'Average error: '
  write (*,'(a,e20.12,a,e20.12)') 'real: ', tmp_r, '    imaginary: ', tmp_i

#if PRINT_ARRAYS
  print *, 'Difference with the original data: '
  do i=0, NC/2-1
     print *, data_complex_copy(i)
  end do
#endif

  print *, '2-------------------------------------------------------------'
  print *, 'Test of the complex array, represented by a real array:'
  do i=0,NC-1
     phase         = 2.0*sll_pi*real(i)/real(NC)
     data_real(i)  = test_func(phase)
  end do

#if PRINT_ARRAYS
  print *, 'Original real data (data_real):'
  do i=0, NC/2-1
     print *, data_real(2*i:2*i+1)
  end do
#endif

  data_real_copy(:) = data_real(:)
  ! interpret data_real as a half-sized sequence of complex numbers:
  do i=0, NC/2-1
     data_complex(i) = transfer(data_real(2*i:),data_complex(0))
  end do
  data_complex_copy(:) = data_complex(:)

#if PRINT_ARRAYS
  print *, 'Original data, as complex, for comparison (data_complex):'
  do i=0, NC/2-1
     print *, data_complex(i)
  end do
  print *, 'transformed complex array, for comparison...'
#endif

  call compute_twiddles(NC/2,twiddles_c_n_2(0:NC/4-1))
  call bit_reverse(NC/4, twiddles_c_n_2(0:NC/4-1))
  call fft_dit_nr_aux( data_complex(0:NC/2-1),   &
                       NC/2,                     &
                       twiddles_c_n_2(0:NC/4-1), &
                       0,                        &
                       FFT_FORWARD )

#if PRINT_ARRAYS
  print *, 'transformed complex array: '
  do i=0,NC/2-1
     print *, data_complex(i)
  end do

 ! call bit_reverse(NC/4, twiddles_c_n_2(0:NC/2-1)) ! natural order now
  print *, 'complex twiddles:'
  do i=0, NC/4-1
     print *, twiddles_c_n_2(i)
  end do
#endif
  print *, 'applying the inverse transform'
  call bit_reverse(NC/4, twiddles_c_n_2(0:NC/4-1)) ! now in natural order
 ! call  fft_dit_nr_aux( data_complex(0:NC/2-1), NC/2, twiddles_c_n_2(0:NC/4-1), 0, FFT_INVERSE )
  call fft_dit_rn_aux( data_complex(0:NC/2-1), NC/2, twiddles_c_n_2(0:NC/4-1),&
                       1, FFT_INVERSE )
  data_complex(:)      = data_complex(:)*1.0/real(NC/2)
  data_complex_copy(:) = data_complex_copy(:) - data_complex(:)
  tmp_c = (0.0,0.0)
  tmp_r = 0.0
  tmp_i = 0.0
  do i=0,size(data_complex_copy)-1
     tmp_r = tmp_r + abs( real(data_complex_copy(i)))
     tmp_i = tmp_i + abs(aimag(data_complex_copy(i)))
  end do
  tmp_r = tmp_r/size(data_complex_copy)
  tmp_i = tmp_i/size(data_complex_copy)
  print *, 'Average error: '
  write (*,'(a,e20.12,a,e20.12)') 'real: ', tmp_r, '    imaginary: ', tmp_i


  print *, 'transformed real array: '
  call compute_twiddles_real_array(NC/2, twiddles_r_n_2(0:NC/2-1))
  call bit_reverse_in_pairs( NC/4, twiddles_r_n_2(0:NC-1))
  call fft_dit_nr_real_array_aux( data_real(0:NC-1),      &
                                  NC/2,                   &
                                  twiddles_r_n_2(0:NC-1), &
                                  0,                      &
                                  FFT_FORWARD )

#if PRINT_ARRAYS
  do i=0, NC/2-1
     write (*,'(f15.8, f15.8)') data_real(2*i), data_real(2*i+1)
  end do
#endif

  call bit_reverse_in_pairs( NC/4, twiddles_r_n_2(0:NC-1) )

#if PRINT_ARRAYS
  print *, 'twiddles represented by a real array:'
  do i=0,NC/4-1
     print *, twiddles_r_n_2(2*i:2*i+1)
  end do
#endif

  print *, 'Entering the inverse real array: '
  call fft_dit_rn_real_array_aux( data_real(0:NC-1),      &
                                  NC/2,                   &
                                  twiddles_r_n_2(0:NC-1), &
                                  1,                      &
                                  FFT_INVERSE )
  data_real(:) = data_real(:)*2.0/real(NC)
  acc = 0.0
  do i=0,size(data_real)-1
     acc = acc + abs(data_real(i) - data_real_copy(i))
  end do

#if PRINT_ARRAYS
  print *, 'Difference with the original data: '
  do i=0, NC-1
     print *, (data_real(i) - data_real_copy(i))
  end do
#endif

  write (*,'(a, e20.12)') 'Average error complex (as real) FFT case = ', acc/NC
  print *, 'Finished test of complex array represented by real data.'

  print *, '3---------------------------------------------------------------'

  print *, 'Test of a real-valued FFT'
  print *, 'For comparison, we run a complex-valued FFT in which the data is the same as the real data we intend to process in the real case.'
  do i=0,NC-1
     phase             = 2.0*sll_pi*real(i)/real(NC)
     val               = test_func(phase)
     data_real(i)      = val
     data_real_copy(i) = val
     data_complex_n(i) = val
  end do

#if PRINT_ARRAYS
  print *, 'data to be transformed:'
  print *, '     real data                            complex data'
  print *, '___________________________________________________________________'
  do i=0, NC-1
     print *, data_real(i), data_complex_n(i)
  end do
#endif

  call compute_twiddles(NC, twiddles_c_n(0:NC/2-1))
  call bit_reverse(NC/2, twiddles_c_n(0:NC/2-1))
  call fft_dit_nr_aux( data_complex_n(0:NC-1), &
                       NC,                     &
                       twiddles_c_n(0:NC/2-1), &
                       0,                      &
                       FFT_FORWARD )
  print *, 'transformed complex array (in natural order): '
  call bit_reverse(NC, data_complex_n(0:NC-1))

#if PRINT_ARRAYS
  do i=0,NC-1
     print *, data_complex_n(i)
  end do
#endif

  ! For the real-valued FFT we need two pairs of twiddles:
  ! - one set to compute complex FFT's of size N/2,
  ! - another set to glue the results of the complex FFT's into the
  !   solution that we want.
  call compute_twiddles_real_array(NC/2, twiddles_r_n_2(0:NC/2-1))
  call bit_reverse_in_pairs( NC/4, twiddles_r_n_2(0:NC-1))
  call compute_twiddles_real_array(NC, twiddles_r_n(0:NC-1))

#if PRINT_ARRAYS
  print *, 'long twiddles:'
  do i=0, size(twiddles_r_n)/2-1
     write (*, '( a, f15.8, a, f15.8, a)') '(', twiddles_r_n(2*i), ', ', &
          twiddles_r_n(2*i+1),')'
  end do
#endif

  print *, 'Entering the real data FFT function'

#if TEST_TIME
  print *, 'timing test:'
  call ETIME(tarray, tresult)
  print *, tresult
  print *, tarray(1)
  print *, tarray(2)
#endif

  call real_data_fft_dit( data_real,      &
                          NC,             &
                          twiddles_r_n_2, &
                          twiddles_r_n,   &
                          FFT_FORWARD )

#if TEST_TIME
  call ETIME(tarray, tresult)
  print *, 'timing test'
  print *, tresult
  print *, tarray(1)
  print *, tarray(2)
#endif



#if PRINT_ARRAYS
  print *, 'After the real FFT:'
  do i=0,NC/2-1
     print *, data_real(2*i:2*i+1)
  end do
#endif

  data_real(:) = data_real(:)*2.0/real(NC)
  print *, 'Proceeding to invert the real transform...'
  call real_data_fft_dit( data_real,      &
                          NC,             &
                          twiddles_r_n_2, &
                          twiddles_r_n,   &
                          FFT_INVERSE )

  acc = 0.0
  do i=0, NC-1
     acc = acc + abs(data_real_copy(i) - data_real(i))
  end do
  print *, 'Averager error:', acc/NC  

#if PRINT_ARRAYS
  do i=0, NC-1
     print *, abs(data_real_copy(i) - data_real(i))
  end do
#endif

  print *, 'REACHED END OF UNIT TEST'

contains
 
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

  
end program fft_test
