
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
#include "sll_splines.h"
#include "sll_constants.h"
  
#define PRINT_SPLINE_COEFFS 0
  use util_constants
  use test_processes_module
  implicit none

  intrinsic :: cos

  sll_int32                              :: ok
 

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
  call test_spline_1d_hrmt( line, test_flag )
  call test_error_flag( test_flag, test_passed, 'test_spline_1d_hrmt')

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
  print *, 'interpolator_tester_1D_hrmt(), cos(x), derivatives test'
  call interpolator_tester_1d_hrmt( &
       mycos, &
       dmycos, &
       interpolate_derivative, &
       0.0_f64, &
       2.0_f64*sll_pi, &
       test_flag )
  call test_error_flag( test_flag, test_passed, 'interpolator_tester_1d_hrmt')


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

  !print *, 'test_passed : ', test_passed
  print *, 'test default slope values, hermite-hermite case'
  call test_2d_spline_hrmt_hrmt_no_slopes( &
       plane3, &
       test_flag )
  call test_error_flag( test_flag, test_passed, &
       'test_2d_spline_hrmt_hrmt_no_slopes')

  call test_2d_cubic_splines_periodic( sincos_prod, &
  !     (/0.0_f64,1.0*sll_pi,0.0_f64,1.0*sll_pi/), 33,33)
  !   (/0.0_f64,2.0_f64*sll_pi,0.0_f64,2.0_f64/), 33,33)
       (/sll_pi,3.0_f64*sll_pi,sll_pi,3.0_f64*sll_pi/), 33,33, test_flag)

  call test_error_flag( test_flag, test_passed, &
       'test_2d_cubic_splines_periodic')

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
       stop
    end if
  end subroutine test_error_flag

  subroutine test_2d_cubic_splines_periodic( func, lims, npts1, npts2, &
    test_flag )

#ifdef STDF95
    sll_real64 :: func
#else
    procedure(fxy) :: func
#endif
    sll_real64, dimension(:), intent(in) :: lims
    sll_int32, intent(in) :: npts1, npts2
    logical, intent(out) :: test_flag
    type(sll_cubic_spline_2d), pointer :: s
    sll_real64                              :: phase_x1, phase_x2
    sll_real64                              :: acc_2D, x1, x2
    sll_real64, allocatable, dimension(:,:) :: data_2d
    sll_real64, allocatable, dimension(:)   :: coordinates_i
    sll_real64, allocatable, dimension(:)   :: coordinates_j
    sll_real64 :: x1min
    sll_real64 :: x1max
    sll_real64 :: x2min
    sll_real64 :: x2max
    sll_int32  :: i,j,err

    test_flag = .true.
    x1min = lims(1)
    x1max = lims(2)
    x2min = lims(3)
    x2max = lims(4)

    SLL_ALLOCATE(data_2d(npts1, npts2), err)        
    SLL_ALLOCATE(coordinates_i(npts1),err)
    SLL_ALLOCATE(coordinates_j(npts2),err)
        
    do i=1, npts1
       coordinates_i(i) = real(i-1,f64)*(x1max-x1min)/real(npts1-1,f64)+x1min
    end do
        
    do j=1, npts2
       coordinates_j(j) = real(j-1,f64)*(x2max-x2min)/real(npts2-1,f64)+x2min
    end do

    do j=1,npts2
       phase_x2 = coordinates_j(j)
       do i=1,npts1
          phase_x1 = coordinates_i(i)
          data_2d(i,j) = func(phase_x1, phase_x2)
       end do
    end do

    s => new_spline_2d( npts1, npts2, x1min, x1max, x2min, x2max, &
         SLL_PERIODIC, SLL_PERIODIC )

    call compute_spline_2d_prdc_prdc( data_2d, s )
    acc_2D = 0.0
    do j=1, npts2
       do i=1, npts1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - interpolate_value_2D(x1,x2,s))
       end do
    end do
    print *, 'Average cumulative error, spline2d, periodic-periodic: ', &
         acc_2D/(NPX1*NPX2)

    if (acc_2D/(NPX1*NPX2)>=1.e-14) then
       print *, 'test periodic cubic 2d spline: error is too big'
       test_flag = .false.
    endif

  end subroutine test_2d_cubic_splines_periodic
end program spline_tester
