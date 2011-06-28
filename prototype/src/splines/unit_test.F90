program spline_tester
#include "sll_working_precision.h"
#include "sll_splines.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use numeric_constants
  implicit none
  
#define NC 32

  sll_int32 :: err
  sll_int32 :: i
  type(sll_spline_1d), pointer :: sp1
  sll_real64, allocatable, dimension(:) :: data
  sll_real64 :: accumulator
  accumulator = 0.0_f64

  print *, 'Spline module unit tester'
  print *, 'allocate data array'
  SLL_ALLOCATE(data(NC+1), err)
  
  print *, 'initialize data array'
  do i=1,NC+1
!data(i) = 1.0
    data(i) = sin((i-1)*sll_pi/real(NC,f64))
  end do
!  print *, 'data: '
!  print *, data(:)
  print *, 'proceed to allocate the spline...'
  sp1 =>  new_spline_1D( data, NC+1, 0.0_f64, sll_pi, PERIODIC_SPLINE )
  
  !  print *, 'Contents of the spline:'
  !  print *, GET_SPLINE_DELTA(sp1)  
  !  print *, GET_SPLINE_XMIN(sp1)
  !  print *, GET_SPLINE_XMAX(sp1)
  print *, 'cumulative errors: '
  do i=1, NC
     accumulator = accumulator + abs(data(i) - interpolate_value(real(i-1,f64)*sll_pi/real(NC,f64), sp1))
     !    print *, accumulator
  end do
  print *, 'average error at the nodes = '
  print *, accumulator/real(NC,f64)
  write (*,'(a,f8.5)')   'original data(0)    = ', data(1)
  write (*,'(a,f20.15)') &
       'interpolated        = ', interpolate_value( 0.0_f64,sp1)
  write (*,'(a,f20.15)')   'original data(NC/4) = ', data(NC/4+1)
  write (*,'(a,f20.15)') &
       'interpolated        = ', interpolate_value( sll_pi/4.0,sp1)
  !  print *, 'spline coefficients: '
  !  print *, sp1%c(:)
  if( accumulator/real(NC,f64) < 1.0e-15 ) then 
     print *, 'PASSED TEST'
  else
     print *, 'FAILED TEST'
  end if
end program spline_tester
