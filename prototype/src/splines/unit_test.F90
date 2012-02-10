
!*****************************************************************************
!
! Selalib      
! Module: unit_test.F90
!
!> @brief 
!> sll_splines unit test
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr), 
!> Eric Sonnendrücker (sonnen@math.unistra.fr)
!> Edwin Chacon-Golcher(chacongolcher@math.unistra.fr)
!                                  
!*****************************************************************************

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
  use util_constants
  use test_processes_module
  implicit none
  
  sll_int32                              :: err, ok
  sll_int32                              :: i, i_test, j_test
  sll_real64, allocatable, dimension(:)  :: knots
  sll_real64, dimension(1:100)           :: spline_vals
  sll_real64                             :: rnd
  sll_real64                             :: reduction
  sll_int32, parameter                   :: nbtest = 12
  logical                                :: test_passed
  print *, 'Test of the 1D spline: '
  ok = 1
  do i_test=1,nbtest     
     call test_process_1d(i_test, ok)     
  enddo
  
  print *, '***************************************************'
  print *, 'Test of the 2D spline: '
  do j_test=1,nbtest
     do i_test=1,nbtest 
        call test_process_2d(i_test, j_test, ok)        
     enddo
  enddo  
  print *, '********************************'


  ! Test with impulse functions
  print *, '***************************************************'
  print *, 'Test of the 1D spline with impulse fonctions: '
  do i_test=13,NP+12     
     call test_process_1d(i_test, ok)     
  enddo

  print *, '***************************************************'
  print *, 'Test of the 2D spline with impulse fonctions: '
  do j_test=13,NPX2+12     
     do i_test=13,NPX1+12
        call test_process_2d(i_test, j_test, ok)
     enddo
  enddo

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

  ! The following tests are currently not included in determining the
  ! value of the OK flag, this should be fixed. The OK flag should be a
  ! logical variable...

  call interpolator_tester_2d( &
       coscos, &
       msincos, &
       interpolate_x1_derivative_2D, &
       test_passed )

  call interpolator_tester_2d( &
       coscos, &
       mcossin, &
       interpolate_x2_derivative_2D, &
       test_passed )

  if (ok==1) then
     print *, ' '
     print *, 'Splines unit test: PASS'
     print *, ' '
  else
     print *, ' '
     print *, 'Splines unit test: FAIL' ! YES BUT WHERE???!!!
     print *, ' '
  endif
  
end program spline_tester
