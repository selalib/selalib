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
  

#define NP 512

  sll_int32 :: err
  sll_int32 :: i, j
  type(sll_spline_1d), pointer :: sp1
  type(sll_spline_1d), pointer :: sp2
  type(sll_spline_2d), pointer :: sp2d
  sll_real64                   :: x
  sll_real64                   :: phase, phase_x1, phase_x2
  sll_real64, allocatable, dimension(:) :: data

  sll_real64, allocatable, dimension(:) :: deriv
  sll_real64, allocatable, dimension(:) :: coordinates
  sll_real64, allocatable, dimension(:) :: out
  sll_real64 :: accumulator1, accumulator2, accumulator3, accumulator4
  sll_real64 :: accumulator5, accumulator6

  ! 2D stuff
  sll_real64  :: acc_2D, x1, x2
  sll_real64, allocatable, dimension(:,:) :: data_2d
  sll_real64, allocatable, dimension(:)   :: coordinates_i
  sll_real64, allocatable, dimension(:)   :: coordinates_j
  ! for specific functions of spline values at a point
  sll_real64, allocatable, dimension(:) :: knots
  sll_real64, dimension(1:100)          :: spline_vals
  sll_real64                            :: rnd
  sll_real64                            :: reduction
#define NPX1  65
#define NPX2  129
#define X1MIN 0.0_f64
#define X1MAX (2.0_f64*sll_pi)
#define X2MIN 0.0_f64
#define X2MAX (4.0_f64*sll_pi)


#define XMIN (-sll_pi)
#define XMAX ( sll_pi)

  accumulator1 = 0.0_f64
  accumulator2 = 0.0_f64
  accumulator3 = 0.0_f64
  accumulator4 = 0.0_f64
  accumulator5 = 0.0_f64
  accumulator6 = 0.0_f64
  print *, 'Spline module unit tester'
  print *, 'allocate data array'
  SLL_ALLOCATE(data(NP), err)
  SLL_ALLOCATE(deriv(NP), err)
  SLL_ALLOCATE(out(NP), err)
  SLL_ALLOCATE(coordinates(NP), err)
  print *, 'initialize data and coordinates array'
  do i=1,NP
     phase          = real(i-1,f64)*(XMAX-XMIN)/real(NP-1,f64) + XMIN
     data(i)        = 2.0*(sin(phase))
     deriv(i)       = 2.0*(cos(phase))
     coordinates(i) = phase
  end do
  !  print *, 'data: '
  !  print *, data(:)
  print *, 'proceed to allocate the spline...'
  sp1 =>  new_spline_1D( NP, XMIN, XMAX, PERIODIC_SPLINE )
  call compute_spline_1D( data, PERIODIC_SPLINE, sp1 )
  sp2 =>  new_spline_1D( NP, XMIN, XMAX, HERMITE_SPLINE, -2.0_f64, -2.0_f64 )
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
  print *, 'periodic case, NP-1 points: '
  print *, 'interpolating individual values from 1 to NP-1:'
  do i=1, NP-1
     x = real(i-1,f64)*(XMAX-XMIN)/real(NP-1,f64)+XMIN
     accumulator1 = accumulator1 + abs(data(i) - interpolate_value(x, sp1))
  end do
  print *, 'checking periodicity:'
  print *, 'difference between values at points 1 and NP: ', &
       abs(data(1) - interpolate_value(XMAX,sp1))
  print *, 'interpolating the whole array:'
  call interpolate_array_values(coordinates, out, NP-1, sp1)
  do i=1, NP-1
     accumulator3 = accumulator3 + abs(data(i) - out(i))
!     write (*, '(a, i8, a, e20.12)') &
!          'point: ', i, ', cumulative err = ',accumulator3
  end do

  print *, 'hermite case, NP points: '
  do i=1, NP
     x = real(i-1,f64)*(XMAX-XMIN)/real(NP-1,f64) + XMIN
     accumulator2 = accumulator2 + abs(data(i) - interpolate_value(x, sp2))
!          sp2%interpolate(real(i-1,f64)*sll_pi/real(NC,f64)))
!     write (*, '(a, i8, a, e20.12)') &
!          'point: ', i, ', cumulative err = ',accumulator2
  end do
  call interpolate_array_values(coordinates, out, NP, sp2)
  do i=1, NP
     accumulator4 = accumulator4 + abs(data(i) - out(i))
!     write (*, '(a, i8, a, e20.12)') &
!          'point: ', i, ', cumulative err = ',accumulator4
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
       'interpolated        = ', interpolate_value( XMIN,sp1)

  write (*,'(a,f20.15)')   'original data((NP-1)/4) = ', data((NP-1)/4)
  write (*,'(a,f20.15)') &
       'interpolated        = ', interpolate_value( (XMAX-XMIN)/4.0+XMIN,sp1)


#if TEST_INTEGRATION
  print *, 'integrating the periodic spline...'
  print *, gauss_legendre_integrate_1D( interpolate_value, sp1, XMIN, XMAX,4)
#endif

#if PRINT_SPLINE_COEFFS
  print *, 'spline coefficients: '
  print *, sp1%coeffs(:)
#endif



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
       'interpolated        = ', interpolate_value( (XMAX-XMIN)/4.0,sp2)
  print *, 'spline coefficients: '
#if PRINT_SPLINE_COEFFS
  print *, sp2%coeffs(:)
#endif

  print *, 'check the slopes, hermite case, first and last points: '
  print *, 'first point: ', (interpolate_value(XMIN+0.0001_f64, sp2) - interpolate_value(XMIN, sp2))/0.0001, '  Expected value= -2.0'
  print *, 'last point: ', (interpolate_value(XMAX, sp2) - interpolate_value(XMAX-0.0001_f64,sp2))/0.0001, '  Expected value= -2.0'



  if( (accumulator2/real(NP,f64) < 1.0e-15) .and. &
      (accumulator4/real(NP,f64) < 1.0e-15) ) then 
     print *, 'PASSED TEST'
  else
     print *, 'FAILED TEST'
  end if

  print *, '---------------------------------------------'
  print *, 'DERIVATIVES TEST'
  do i=1, NP
     phase          = real(i-1,f64)*(XMAX-XMIN)/real(NP-1,f64)+XMIN
     accumulator5 = accumulator5 + abs(deriv(i) - &
          interpolate_derivative(phase, sp1))
     accumulator6 = accumulator6 + abs(deriv(i) - &
          interpolate_derivative(phase, sp2))
 !    write (*,'(a, f20.15, a, f20.15, a, f20.15, a, f20.15)') 'phase = ', phase, ', derivative value = ', deriv(i), ', computed value (periodic)= ', interpolate_derivative(phase, sp1), ', computed value (hermite) = ', interpolate_derivative(phase,sp2)
  end do
  print *, 'average error at the nodes (single values, periodic) = '
  print *, accumulator5/real(NP,f64)
  print *, 'average error at the nodes (single values, hermite) = '
  print *, accumulator6/real(NP,f64)
  if( &
     (accumulator5/real(NP,f64) .le. 1.0e-6) .and. &
     (accumulator6/real(NP/f64) .le. 1.0e-6) ) then
     print *, 'PASSED DERIVATIVES TEST'
  else
     print *, 'FAILED DERIVATIVES TEST' 
  end if
  call delete_spline_1D(sp1)
  call delete(sp2)
  print *, '***************************************************'
  print *, 'Test of the 2D spline: '
  SLL_ALLOCATE(data_2d(NPX1, NPX2), err)
  print *, 'Filling data:'

  SLL_ALLOCATE(coordinates_i(NPX1),err)
  SLL_ALLOCATE(coordinates_j(NPX2),err)

  do j=1,NPX2
     do i=1,NPX1
        phase_x1 = real(i-1,f64)*X1MAX/real(NPX1-1,f64)
        phase_x2 = real(j-1,f64)*X2MAX/real(NPX2-1,f64)
        coordinates_i(i) = phase_x1
        coordinates_j(j) = phase_x2
        data_2d(i,j) = sin(phase_x1)*sin(phase_x2)
     end do
  end do
  print *, 'Allocating 2D spline...'

  sp2d => new_spline_2D(NPX1, NPX2, X1MIN, X1MAX, X2MIN, X2MAX, PERIODIC_SPLINE, PERIODIC_SPLINE)

  print *, 'Computing the 2D spline...'
  call compute_spline_2D_prdc_prdc( data_2d, sp2d )
  print *, 'Completed computing the 2d spline.'

  acc_2D = 0.0
  do j=1, NPX2
     do i=1, NPX1
        x1 = coordinates_i(i)
        x2 = coordinates_j(j)
        acc_2D = acc_2D + &
             abs(data_2d(i,j) - interpolate_value_2D(x1,x2,sp2d))
     end do
  end do
  print *, 'Cumulative error, spline2d, periodic-periodic: ', acc_2D
#if 0
  print *, 'Deleting the 2D spline...'
  call delete(sp2d)
#endif

  print *, '********************************'
  print *, 'spline values at points:'
#define NUM_KNOTS 10
#define SPLINE_DEGREE 3
  SLL_ALLOCATE(knots(NUM_KNOTS),err)
  knots(1) = 0.0
  do i=2,NUM_KNOTS
     call random_number(rnd)
     knots(i) = knots(i-1) + rnd 
  end do
  print *, 'knots: '
  print *, knots(:)
  spline_vals(1:SPLINE_DEGREE+1) = b_splines_at_x(knots,SPLINE_DEGREE,5,3.3_f64)
  print *, spline_vals(1:SPLINE_DEGREE+1)
  reduction = 0.0

  do i=1,SPLINE_DEGREE+1
     reduction = reduction + spline_vals(i)
  end do
  print *, 'sum of spline values = ', reduction
#undef NUM_KNOTS
#undef SPLINE_DEGREE
  print *, 'END TEST'
end program spline_tester
