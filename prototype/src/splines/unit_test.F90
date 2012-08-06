
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
!> Eric SONNENDRÜCKER (sonnen@math.unistra.fr)
!> Edwin CHACON-GOLCHER (chacongolcher@math.unistra.fr)
!                                  
!*****************************************************************************

program spline_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  
#define PRINT_SPLINE_COEFFS 0
  use sll_splines
  use numeric_constants
  use util_constants
  use test_processes_module
  implicit none

  intrinsic :: cos

  sll_int32                              :: err, ok
  sll_int32                              :: i, i_test, j_test
 
  sll_real64, dimension(1:100)           :: spline_vals
  sll_real64                             :: rnd

  sll_int32, parameter                   :: nbtest = 12
  logical                                :: test_passed
  logical                                :: test_flag
  print *, 'Test of the 1D spline: '
  ok = 1
  test_passed = .true.

#if 0
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
  print *, 'Test of the 1D spline with impulse functions: '
  do i_test=13,NP+12     
     call test_process_1d(i_test, ok)     
  enddo

  print *, '***************************************************'
  print *, 'Test of the 2D spline with impulse functions: '
  do j_test=13,NPX2+12     
     do i_test=13,NPX1+12
        call test_process_2d(i_test, j_test, ok)
     enddo
  enddo

#endif

  ! The following tests are currently not included in determining the
  ! value of the OK flag, this should be fixed. The OK flag should be a
  ! logical variable...
  print *, 'test_spline_1d_hrmt, linear function case: '
  call test_spline_1d_hrmt( line, 1.0_f64, 1.0_f64, test_flag )
  call test_error_flag( test_flag, test_passed, 'test_spline_1d_hrmt')

  call test_spline_1d_hrmt_no_slopes( line, test_flag )
  call test_error_flag( test_flag, test_passed,'test_spline_1d_hrmt_no_slopes')

  print *, ' '
  print *, 'interpolator_tester_1D_prdc(), cos(x), normal values'
  call interpolator_tester_1d_prdc( &
       mycos, &
       mycos, &
       interpolate_value, &
       0.0_f64, &
       2.0_f64*sll_pi, &
       33, &
       test_flag )
  call test_error_flag( test_flag, test_passed, 'interpolator_tester_1d_prdc')

  print *, ' '
  print *, 'interpolator_tester_1D_prdc(), cos(x), derivatives test'
  call interpolator_tester_1d_prdc( &
       mycos, &
       dmycos, &
       interpolate_derivative, &
       0.0_f64, &
       2.0_f64*sll_pi, &
       33, &
       test_flag, &
       6.0e-6_f64)
  call test_error_flag( test_flag, test_passed, 'interpolator_tester_1d_prdc')

  print *, ' '
  print *, 'interpolator_tester_2d(), cos(x)*cos(y) case, deriv in X1: '
  call interpolator_tester_2d_prdc_prdc( &
       coscos, &
       msincos, &
       interpolate_x1_derivative_2D, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'interpolator_tester_2d_prdc_prdc')

  print *, ' '
  print *, 'interpolator_tester_2d(), cos(x)*cos(y) case, deriv in X2: '
  call interpolator_tester_2d_prdc_prdc( &
       coscos, &
       mcossin, &
       interpolate_x2_derivative_2D, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'interpolator_tester_2d_prdc_prdc')

  print *, '----------------------------------------------------'
  print *,  'Test polar transformation case: '
  print *, '----------------------------------------------------'
  print *, ' '
  print *, 'hrmt_prdc on polar_x:'
  call test_2d_spline_hrmt_prdc( &
       polar_x, &
       deriv1_polar_x, &
       deriv1_polar_x, &
       test_flag )
  call test_error_flag( test_flag, test_passed, 'test_2d_spline_hrmt_prdc')

  print *, ' '
  print *, 'hrmt_prdc on polar_x, test default slopes:'
  call test_2d_spline_hrmt_prdc_no_slopes( &
       polar_x, &
       deriv1_polar_x, &
       deriv1_polar_x, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'test_2d_spline_hrmt_prdc_no_slopes')

  print *, ' '
  print *, 'hrmt_prdc on polar_y, deriv1:'
  call test_2d_spline_hrmt_prdc( &
       polar_y, &
       deriv1_polar_y, &
       deriv1_polar_y, &
       test_flag )
  call test_error_flag( test_flag, test_passed, 'test_2d_spline_hrmt_prdc')

  print *, ' '
  print *, 'hrmt_prdc on polar_y, deriv2:'
  call test_2d_spline_hrmt_prdc( &
       polar_y, &
       deriv2_polar_y, &
       deriv2_polar_y, &
       test_flag )
  call test_error_flag( test_flag, test_passed, 'test_2d_spline_hrmt_prdc')

  print *, ' '
  print *, 'interpolator tester, hrmt_prdc, on polar_x, deriv1:'
  call interpolator_tester_2d_hrmt_prdc( &
       polar_x, &
       deriv1_polar_x, &
       interpolate_x1_derivative_2D, &
       deriv1_polar_x, &
       deriv1_polar_x, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'interpolator_tester_2d_hrmt_prdc')

  print *, ' '
  print *, 'interpolator tester, hrmt_prdc, on polar_x, deriv2:'
  call interpolator_tester_2d_hrmt_prdc( &
       polar_x, &
       deriv2_polar_x, &
       interpolate_x2_derivative_2D, &
       deriv1_polar_x, &
       deriv1_polar_x, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'interpolator_tester_2d_hrmt_prdc')

  print *, ' '
  print *, 'interpolator tester, hrmt_prdc, on polar_y, deriv1:'
  call interpolator_tester_2d_hrmt_prdc( &
       polar_y, &
       deriv1_polar_y, &
       interpolate_x1_derivative_2D, &
       deriv1_polar_y, &
       deriv1_polar_y, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'interpolator_tester_2d_hrmt_prdc')

  print *, ' '
  print *, 'interpolator tester, hrmt_prdc, on polar_y, deriv2:'
  call interpolator_tester_2d_hrmt_prdc( &
       polar_y, &
       deriv2_polar_y, &
       interpolate_x2_derivative_2D, &
       deriv1_polar_y, &
       deriv1_polar_y, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'interpolator_tester_2d_hrmt_prdc')


  print *, 'plane test case'
  call test_2d_spline_hrmt_prdc( plane, plane_deriv, plane_deriv, test_flag )
  call test_error_flag( test_flag, test_passed, &
       'test_2d_spline_hrmt_prdc')

  call test_2d_spline_prdc_hrmt( plane2,plane2_deriv,plane2_deriv,test_flag )
  call test_error_flag( test_flag, test_passed, &
       'test_2d_spline_prdc_hrmt')

  print *, ' '
  print *, 'testing the default, computed slope values'
  call test_2d_spline_prdc_hrmt_no_slopes( &
       plane2, &
       plane2_deriv, &
       plane2_deriv, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'test_2d_spline_prdc_hrmt_no_slopes')


  call test_2d_spline_hrmt_hrmt( &
       plane3, &
       plane3_deriv_x, &
       plane3_deriv_x, &
       plane3_deriv_y, &
       plane3_deriv_y, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'test_2d_spline_hrmt_hrmt')

  print *, ' '
  print *, 'test default slope values, hermite-hermite case'
  call test_2d_spline_hrmt_hrmt_no_slopes( &
       plane3, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'test_2d_spline_hrmt_hrmt_no_slopes')



  if (test_passed .eqv. .true.) then
     print *, ' '
     print *, 'Splines unit test: PASSED'
     print *, ' '
  else
     print *, ' '
     print *, 'Splines unit test: FAILED'
     print *, ' '
  endif
  
contains

  subroutine test_error_flag( individual_flag, general_flag, message )
    logical, intent(in) :: individual_flag 
    logical, intent(inout) :: general_flag
    character(len=*) :: message
    ! print *, individual_flag
    general_flag = general_flag .and. individual_flag
    if( individual_flag .eqv. .false. ) then
       print *, 'FAILURE IN FUNCTION: ', message
    end if
  end subroutine test_error_flag

end program spline_tester
