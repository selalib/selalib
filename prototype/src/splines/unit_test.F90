program spline_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

#define TEST_INTEGRATION 0
#if TEST_INTEGRATION
  use gauss_legendre_integration
#endif
#define PRINT_SPLINE_COEFFS 0
  use sll_splines
  use numeric_constants
  implicit none
  

#define NP 5000

  sll_int32 :: err
  sll_int32 :: i
  type(sll_spline_1d), pointer :: sp1
  type(sll_spline_1d), pointer :: sp2
  sll_real64                   :: phase
  sll_real64, allocatable, dimension(:) :: data
  sll_real64, allocatable, dimension(:) :: coordinates
  sll_real64, allocatable, dimension(:) :: out
  sll_real64 :: accumulator1, accumulator2, accumulator3, accumulator4

  accumulator1 = 0.0_f64
  accumulator2 = 0.0_f64
  accumulator3 = 0.0_f64
  accumulator4 = 0.0_f64
  print *, 'Spline module unit tester'
  print *, 'allocate data array'
  SLL_ALLOCATE(data(NP), err)
  SLL_ALLOCATE(out(NP), err)
  SLL_ALLOCATE(coordinates(NP), err)
  print *, 'initialize data and coordinates array'
  do i=1,NP
     phase          = real(i-1,f64)*sll_pi/real(NP-1,f64)
     data(i)        = sin(phase)
     coordinates(i) = phase
  end do
  !  print *, 'data: '
  !  print *, data(:)
  print *, 'proceed to allocate the spline...'
  sp1 =>  new_spline_1D( NP, 0.0_f64, sll_pi, PERIODIC_SPLINE )
  call compute_spline_1D( data, PERIODIC_SPLINE, sp1 )
  sp2 =>  new_spline_1D( NP, 0.0_f64, sll_pi, HERMITE_SPLINE )
  call compute_spline_1D( data, HERMITE_SPLINE, sp2 )

  print *, 'Contents of the spline 1:'
  print *, sp1%xmin
  print *, sp1%xmax
  print *, sp1%delta
  print *, sp1%rdelta
  print *, sp1%bc_type
#if PRINT_SPLINE_COEFFS
  print *, sp1%coeffs(:)
#endif
  print *, 'Contents of the spline 2:'
  print *, sp2%xmin
  print *, sp2%xmax
  print *, sp2%delta
  print *, sp2%rdelta
  print *, sp2%bc_type
#if PRINT_SPLINE_COEFFS
  print *, sp2%coeffs(:)
#endif
  print *, 'cumulative errors: '
  print *, 'periodic case, NP points: '
  print *, 'interpolating individual values:'
  do i=1, NP
     accumulator1 = accumulator1 + abs(data(i) - &
          interpolate_value(real(i-1,f64)*sll_pi/real(NP-1,f64), sp1))
!         sp1%interpolate(real(i-1,f64)*sll_pi/real(NC,f64)))
     write (*, '(a, i8, a, e20.12)') &
          'point: ', i, ', cumulative err = ',accumulator1
  end do
  print *, 'interpolating the whole array:'
  call interpolate_array_values(coordinates, out, NP, sp1)
  do i=1, NP
     accumulator3 = accumulator3 + abs(data(i) - out(i))
     write (*, '(a, i8, a, e20.12)') &
          'point: ', i, ', cumulative err = ',accumulator3
  end do

  print *, 'hermite case, NP points: '
  do i=1, NP
     accumulator2 = accumulator2 + abs(data(i) - &
          interpolate_value(real(i-1,f64)*sll_pi/real(NP-1,f64), sp2))
!          sp2%interpolate(real(i-1,f64)*sll_pi/real(NC,f64)))
     write (*, '(a, i8, a, e20.12)') &
          'point: ', i, ', cumulative err = ',accumulator2
  end do
  call interpolate_array_values(coordinates, out, NP, sp2)
  do i=1, NP
     accumulator4 = accumulator4 + abs(data(i) - out(i))
     write (*, '(a, i8, a, e20.12)') &
          'point: ', i, ', cumulative err = ',accumulator4
  end do
  print *, '----------------------------------------------------'
  print *, 'RESULTS: '
  print *, 'Periodic case: '
  print *, 'average error at the nodes (single values) = '
  print *, accumulator1/real(NP,f64)
  print *, 'average error at the nodes (whole array) = '
  print *, accumulator3/real(NP,f64)
  write (*,'(a,f8.5)')   'original data(0)    = ', data(1)
  write (*,'(a,f20.15)') &
       'interpolated        = ', interpolate_value( 0.0_f64,sp1)
  write (*,'(a,f20.15)')   'original data(NC/4) = ', data((NP-1)/4+1)
  write (*,'(a,f20.15)') &
       'interpolated        = ', interpolate_value( sll_pi/4.0,sp1)

#if TEST_INTEGRATION
  print *, 'integrating the periodic spline...'
  print *, gauss_legendre_integrate_1D( interpolate_value, sp1, 0.0_f64, sll_pi,4)
#endif

#if PRINT_SPLINE_COEFFS
  print *, 'spline coefficients: '
  print *, sp1%coeffs(:)
#endif
  call delete_spline_1D(sp1)
  if( (accumulator1/real(NP,f64) < 1.0e-15) .and. &
      (accumulator3/real(NP,f64) < 1.0e-15) ) then 
     print *, 'PASSED TEST'
  else
     print *, 'FAILED TEST'
  end if
  print *, '**************************** '
  print *, 'Hermite case: '
  print *, 'average error at the nodes = '
  print *, accumulator2/real(NP,f64)
  write (*,'(a,f8.5)')   'original data(0)    = ', data(1)
  write (*,'(a,f20.15)') &
       'interpolated        = ', interpolate_value( 0.0_f64,sp2)
  write (*,'(a,f20.15)')   'original data((NP-1)/4) = ', data((NP-1)/4+1)
  write (*,'(a,f20.15)') &
       'interpolated        = ', interpolate_value( sll_pi/4.0,sp2)
  print *, 'spline coefficients: '
#if PRINT_SPLINE_COEFFS
  print *, sp2%coeffs(:)
#endif
  call delete(sp2)
  if( (accumulator2/real(NP,f64) < 1.0e-15) .and. &
      (accumulator4/real(NP,f64) < 1.0e-15) ) then 
     print *, 'PASSED TEST'
  else
     print *, 'FAILED TEST'
  end if


end program spline_tester
