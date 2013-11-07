program unit_test
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_interpolators.h"
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_gnuplot
  implicit none

#define NPTS1 64
#define NPTS2 64 
#define SPL_DEG 9
#define X1MIN 0.0_f64
#define X1MAX 1.0_f64
#define X2MIN 0.0_f64
#define X2MAX 1.0_f64
#define TOLERANCE_NODE 1.0E-7_f64
#define TOLERANCE_DER  3.0e-5_f64

  type(arb_deg_2d_interpolator) :: ad2d
  sll_real64, dimension(:,:), allocatable    :: x
  sll_real64, dimension(:), allocatable      :: eta1_pos
  sll_real64, dimension(:), allocatable      :: eta2_pos
  sll_real64, dimension(:,:), allocatable    :: reference
  sll_real64, dimension(:,:), allocatable    :: calculated
  sll_real64, dimension(:,:), allocatable    :: difference
!!$  sll_real64, dimension(:), allocatable      :: x1_eta1_min
!!$  sll_real64, dimension(:), allocatable      :: x1_eta1_max
  sll_int32 :: ierr
  sll_int32  :: i, j
  sll_real64 :: eta1, eta2, h1, h2
  sll_real64 :: acc, acc1, acc2, acc3, node_val, ref, deriv1_val 
  sll_real64 :: deriv2_val
  sll_real64 :: acc_der1, acc1_der1, acc2_der1, acc3_der1
  sll_real64 :: acc_der2, acc1_der2, acc2_der2, acc3_der2
  sll_real64 :: normL2_0, normL2_1, normL2_2, normL2_3
  sll_real64 :: normH1_0, normH1_1, normH1_2, normH1_3
  logical :: result

  
  result = .true.
  print *,  'filling out discrete arrays for x1 '
  h1 = (X1MAX-X1MIN)/real(NPTS1-1,f64)
  h2 = (X2MAX-X2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  allocate(x(NPTS1,NPTS2))
  allocate(reference(NPTS1,NPTS2))
  allocate(calculated(NPTS1,NPTS2))
  allocate(difference(NPTS1,NPTS2))
  allocate(eta1_pos(NPTS1))
  allocate(eta2_pos(NPTS2))
!  allocate(x1_eta1_min(NPTS2))
 ! allocate(x1_eta1_max(NPTS2))
  print *, '***********************************************************'
  print *, '              periodic-periodic case'
  print *, '***********************************************************'
  
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1               = X1MIN + real(i,f64)*h1
        eta2               = X2MIN + real(j,f64)*h2
        eta1_pos(i+1)      = eta1
        eta2_pos(j+1)      = eta2
        x(i+1,j+1)         = cos(2.0_f64*sll_pi*eta1)!cos(2.0_f64*sll_pi*eta2) *cos(2.0_f64*sll_pi*eta1)
        reference(i+1,j+1) = cos(2.0_f64*sll_pi*eta1)!cos(2.0_f64*sll_pi*eta2)*cos(2.0_f64*sll_pi*eta1)
     end do
  end do
  
  call sll_gnuplot_field_2d( &
       X1MIN, &
       X1MAX, &
       NPTS1, &
       X2MIN,&
       X2MAX,&
       NPTS2, &
       reference, &
       'reference_interp_arb_deg', &
       0, &
       ierr)
  
  
  ! Test the 2D transformation:
  
  call ad2d%initialize( &
       NPTS1, &
       NPTS2, &
       X1MIN, &
       X1MAX, &
       X2MIN, &
       X2MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SPL_DEG, &
       SPL_DEG )
  

  call ad2d%compute_interpolants( &
       x(1:NPTS1,1:NPTS2),&
       eta1_pos(1:NPTS1),&
       NPTS1,&
       eta2_pos(1:NPTS2),&
       NPTS2)
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  acc      = 0.0_f64
  acc_der1 = 0.0_f64
  acc_der2 = 0.0_f64
  normL2_0 = 0.0_f64
  normH1_0 = 0.0_f64

  do j=0,NPTS2-2
     do i=0,NPTS1-2
        eta1       = X1MIN + real(i,f64)*h1
        eta2       = X2MIN + real(j,f64)*h2
        node_val   = ad2d%interpolate_value(eta1,eta2)
        ref        = cos(2.0_f64*sll_pi*eta1)!cos(2.0_f64*sll_pi*eta2)*cos(2.0_f64*sll_pi*eta1)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        !print*, eta1,eta2,node_val,ref,ref-node_val
        
        acc        = acc + abs(node_val-ref)

        normL2_0 = normL2_0 + (node_val-ref)**2 *h1*h2
        deriv1_val = ad2d%interpolate_derivative_eta1(eta1,eta2)
        ref = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)!cos(2.0_f64*sll_pi*eta2)*sin(2.0_f64*sll_pi*eta1)
        acc_der1 = acc_der1 + abs(deriv1_val-ref)
        !
        !print*,'derive=', ref,deriv1_val,ref-deriv1_val
        normH1_0 = normH1_0 + (deriv1_val-ref)**2 *h1*h2
        deriv2_val = ad2d%interpolate_derivative_eta2(eta1,eta2)
        ref  = 0.0_f64!-2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta2)*cos(2.0_f64*sll_pi*eta1)
        acc_der2 = acc_der2 + abs(deriv2_val-ref)
        !print*, ref,deriv2_val
        normH1_0 = normH1_0 + (deriv2_val-ref)**2 *h1*h2
        !print*, deriv1_val,deriv2_val
     end do
  end do

  call sll_gnuplot_field_2d(X1MIN, X1MAX, NPTS1, X2MIN, X2MAX, NPTS2, &
       calculated, 'calculated_interp_arb_deg', 0, ierr)

  call sll_gnuplot_field_2d(X1MIN, X1MAX, NPTS1, X2MIN, X2MAX, NPTS2, &
       difference, 'difference_interp_arb_deg', 0, ierr)

  call sll_gnuplot_field_2d(X1MIN, X1MAX, NPTS1, X2MIN ,X2MAX, NPTS2, &
       ad2d%coeff_splines, 'coefficients_interp_arb_deg', 0, ierr)
  
  
  print *, '***********************************************************'
  print *, '              periodic-dirichlet case'
  print *, '***********************************************************'

  call delete(ad2d)

  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1               = X1MIN + real(i,f64)*h1
        eta2               = X2MIN + real(j,f64)*h2
        eta1_pos(i+1)      = eta1
        eta2_pos(j+1)      = eta2
        x(i+1,j+1)         = sin(2.0_f64*sll_pi*eta2) *sin(2.0_f64*sll_pi*eta1)
        reference(i+1,j+1) = sin(2.0_f64*sll_pi*eta2)*sin(2.0_f64*sll_pi*eta1)
     end do
  end do
  
  
  call ad2d%initialize( &
       NPTS1, &
       NPTS2, &
       X1MIN, &
       X1MAX, &
       X2MIN, &
       X2MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPL_DEG, &
       SPL_DEG )

  call ad2d%compute_interpolants( &
       x(1:NPTS1,1:NPTS2),&
       eta1_pos(1:NPTS1),&
       NPTS1,&
       eta2_pos(1:NPTS2),&
       NPTS2)

  
  print *, 'Compare the values of the transformation at the nodes: '

  acc1 = 0.0_f64
  acc1_der1 = 0.0_f64
  acc1_der2 = 0.0_f64
  normL2_1  = 0.0_F64
  normH1_1  = 0.0_f64
  do j=0,NPTS2-1
     do i=0,NPTS1-2
        eta1       = X1MIN + real(i,f64)*h1
        eta2       = X2MIN + real(j,f64)*h2
        node_val   = ad2d%interpolate_value(eta1,eta2)
        ref        = sin(2.0_f64*sll_pi*eta2)*sin(2.0_f64*sll_pi*eta1)
        normL2_1 = normL2_1 + (node_val-ref)**2 *h1*h2
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        !print*, ref,node_val,node_val-ref
        !print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
         !    'theoretical = ', ref,'difference=',ref-node_val
        acc1        = acc1 + abs(node_val-ref)
        
        deriv1_val = ad2d%interpolate_derivative_eta1(eta1,eta2)   
        ref = 2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta2)*cos(2.0_f64*sll_pi*eta1)
        acc1_der1 = acc1_der1 + abs(deriv1_val-ref)
        !print*, ref,deriv1_val
        normH1_1 = normH1_1 + (deriv1_val-ref)**2 *h1*h2
        deriv2_val = ad2d%interpolate_derivative_eta2(eta1,eta2)
        ref  = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta2)*sin(2.0_f64*sll_pi*eta1)
        acc1_der2 = acc1_der2 + abs(deriv2_val-ref)
        !print*, ref,deriv2_val
        normH1_1 = normH1_1 + (deriv2_val-ref)**2 *h1*h2
     
     end do
  end do
  
  
  
  print *, '***********************************************************'
  print *, '              dirichlet-periodic case'
  print *, '***********************************************************'
  
  call delete(ad2d)

  !reinitialize data
  ! assumes eta mins are 0
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1               = X1MIN + real(i,f64)*h1
        eta2               = X2MIN + real(j,f64)*h2
        eta1_pos(i+1)      = eta1
        eta2_pos(j+1)      = eta2
        x(i+1,j+1)         = sin(2.0_f64*sll_pi*eta1)*cos(2.0_f64*sll_pi*eta2)
        reference(i+1,j+1) = sin(2.0_f64*sll_pi*eta1)
     end do
  end do

  call ad2d%initialize( &
       NPTS1, &
       NPTS2, &
       X1MIN, &
       X1MAX, &
       X2MIN, &
       X2MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SPL_DEG, &
       SPL_DEG )

  call ad2d%compute_interpolants( &
       x(1:NPTS1,1:NPTS2),&
       eta1_pos(1:NPTS1),&
       NPTS1,&
       eta2_pos(1:NPTS2),&
       NPTS2)
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  acc2 = 0.0_f64
  acc2_der1 = 0.0_f64
  acc2_der2 = 0.0_f64
  normL2_2  = 0.0_F64
  normH1_2  = 0.0_f64
  do j=0,NPTS2-2
     do i=0,NPTS1-1
        eta1       = X1MIN + real(i,f64)*h1
        eta2       = X2MIN + real(j,f64)*h2
        !print*, "hehe"
        node_val   = ad2d%interpolate_value(eta1,eta2)
        !print*, "hehe"
        ref                 = sin(2.0_f64*sll_pi*eta1)*cos(2.0_f64*sll_pi*eta2)
        normL2_2 = normL2_2 + (node_val-ref)**2 *h1*h2
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        !print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
         !    'theoretical = ', ref
        acc2        = acc2 + abs(node_val-ref)
        deriv1_val = ad2d%interpolate_derivative_eta1(eta1,eta2)
        
        ref = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)*cos(2.0_f64*sll_pi*eta2)
        acc2_der1 = acc2_der1 + abs(deriv1_val-ref)
        !print*, ref,deriv1_val
        normH1_2 = normH1_2 + (deriv1_val-ref)**2 *h1*h2
        deriv2_val = ad2d%interpolate_derivative_eta2(eta1,eta2)
        ref  = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)*sin(2.0_f64*sll_pi*eta2)
        acc2_der2 = acc2_der2 + abs(deriv2_val-ref)
        !print*, ref,deriv2_val
        normH1_2 = normH1_2 + (deriv2_val-ref)**2 *h1*h2
     end do
  end do
  
  print *, '***********************************************************'
  print *, '              dirichlet-dirichlet case'
  print *, '***********************************************************'
  
  call delete(ad2d)

  !reinitialize data

  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1               = X1MIN + real(i,f64)*h1
        eta2               = X2MIN + real(j,f64)*h2
        eta1_pos(i+1)      = eta1
        eta2_pos(j+1)      = eta2
        x(i+1,j+1)         = sin(2.0_f64*sll_pi*eta1)*sin(2.0_f64*sll_pi*eta2) 
        reference(i+1,j+1) = sin(2.0_f64*sll_pi*eta1)*sin(2.0_f64*sll_pi*eta2)
     end do
  end do
  
  call ad2d%initialize( &
       NPTS1, &
       NPTS2, &
       X1MIN, &
       X1MAX, &
       X2MIN, &
       X2MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPL_DEG, &
       SPL_DEG )
  !print*, 'ret', size(x,1), size(x,2)
  !print*, 'x', x(1:NPTS1,1:NPTS2)
  call ad2d%compute_interpolants( &
       x,&
       eta1_pos,&
       NPTS1,&
       eta2_pos,&
       NPTS2)

  
  !node_val   = ad2d%interpolate_value(0.0_f64,0.0_f64)
  print *, 'Compare the values of the transformation at the nodes: '

  acc3 = 0.0_f64
  acc3_der1 = 0.0_f64
  acc3_der2 = 0.0_f64
  normL2_3  = 0.0_F64
  normH1_3  = 0.0_f64
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1       = X1MIN + real(i,f64)*h1
        eta2       = X2MIN + real(j,f64)*h2
        !print*, "hehe"
        node_val   = ad2d%interpolate_value(eta1,eta2)
        !print*, "hehe"
        ref                 = reference(i+1,j+1)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        normL2_3 = normL2_3 + (node_val-ref)**2 *h1*h2
        !print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
         !    'theoretical = ', ref
        acc3        = acc3 + abs(node_val-ref)

        deriv1_val = ad2d%interpolate_derivative_eta1(eta1,eta2)        
        ref = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)*sin(2.0_f64*sll_pi*eta2)
        acc3_der1 = acc3_der1 + abs(deriv1_val-ref)
        !print*, ref,deriv1_val,abs(deriv1_val-ref)
        normH1_3 = normH1_3 + (deriv1_val-ref)**2 *h1*h2
        deriv2_val = ad2d%interpolate_derivative_eta2(eta1,eta2)
        ref  = 2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)*cos(2.0_f64*sll_pi*eta2)
        acc3_der2 = acc3_der2 + abs(deriv2_val-ref)
        !print*, ref,deriv2_val
        normH1_3 = normH1_3 + (deriv2_val-ref)**2 *h1*h2
     end do
  end do


  
  print*, '--------------------------------------------'
  print*, ' Average error in nodes'
  print*, '--------------------------------------------'
  print *, 'Average error in nodes (dirichlet-dirichlet) = ', acc3/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc3/(NPTS1*NPTS2), TOLERANCE_NODE,&
       result)
  print *, 'Average error in nodes (dirichlet-periodic) = ', acc2/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc2/(NPTS1*NPTS2), TOLERANCE_NODE, &
       result)
  print *, 'Average error in nodes (periodic-dirichlet) = ', acc1/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc1/(NPTS1*NPTS2), TOLERANCE_NODE, &
       result)
  print *, 'Average error in nodes (periodic-periodic) = ', acc/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc/(NPTS1*NPTS2), TOLERANCE_NODE, &
       result)

  print*, '--------------------------------------------'
  print*, ' Average error in nodes first derivative eta1'
  print*, '--------------------------------------------'
  print *,'Average error in nodes first derivative eta1(dirichlet-dirichlet)=',&
       acc3_der1/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc3_der1/(NPTS1*NPTS2),TOLERANCE_DER,&
       result)
  print *,'Average error in nodes first derivative eta1(dirichlet-periodic)=',&
       acc2_der1/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc2_der1/(NPTS1*NPTS2),TOLERANCE_DER,&
       result)
  print *,'Average error in nodes first derivative eta1(periodic-dirichlet)=',&
       acc1_der1/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc1_der1/(NPTS1*NPTS2),TOLERANCE_DER,&
       result)
  print *,'Average error in nodes first derivative eta1(periodic-periodic)=',&
       acc_der1/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc_der1/(NPTS1*NPTS2),TOLERANCE_DER,&
       result)

  print*, '--------------------------------------------'
  print*, ' Average error in nodes first derivative eta2'
  print*, '--------------------------------------------'
  print *,'Average error in nodes first derivative eta2(dirichlet-dirichlet)=',&
       acc3_der2/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc3_der2/(NPTS1*NPTS2),TOLERANCE_DER,&
       result)
  print *,'Average error in nodes first derivative eta2(dirichlet-periodic)=',&
       acc2_der2/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc2_der2/(NPTS1*NPTS2),TOLERANCE_DER,&
       result)
  print *,'Average error in nodes first derivative eta2(periodic-dirichlet)=',&
       acc1_der2/(NPTS1*NPTS2)
  call test_value_for_acceptable_error( acc1_der2/(NPTS1*NPTS2),TOLERANCE_DER,&
       result)
  print *,'Average error in nodes first derivative eta2(periodic-periodic)=',&
       acc_der2/(NPTS1*NPTS2)

  print*, '--------------------------------------------'
  print*, ' Error norm L2'
  print*, '--------------------------------------------'
  print *,'Error norm L2 (dirichlet-dirichlet)=',sqrt(normL2_3), h1**(SPL_DEG)*(2.0_f64*sll_pi)
  print *,'Error norm L2 (dirichlet-periodic)=', sqrt(normL2_2), h1**(SPL_DEG)*(2.0_f64*sll_pi)
  print *,'Error norm L2 (periodic-dirichlet)=', sqrt(normL2_1), h1**(SPL_DEG)*(2.0_f64*sll_pi)
  print *,'Error norm L2 (periodic-periodic)=',  sqrt(normL2_0), h1**(SPL_DEG)*(2.0_f64*sll_pi)

  print*, '--------------------------------------------'
  print*, ' Error norm H1'
  print*, '--------------------------------------------'
  print *,'Error norm H1 (dirichlet-dirichlet)=',sqrt(normH1_3), h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2
  print *,'Error norm H1 (dirichlet-periodic)=', sqrt(normH1_2), h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2
  print *,'Error norm H1 (periodic-dirichlet)=', sqrt(normH1_1), h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2
  print *,'Error norm H1 (periodic-periodic)=',  sqrt(normH1_0), h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2

  if (  ( sqrt(normL2_0) <= h1**(SPL_DEG)*(2.0_f64*sll_pi))   .AND. &
        ( sqrt(normL2_1) <= h1**(SPL_DEG)*(2.0_f64*sll_pi))   .AND. &
        ( sqrt(normL2_2) <= h1**(SPL_DEG)*(2.0_f64*sll_pi))   .AND. &
        ( sqrt(normL2_3) <= h1**(SPL_DEG)*(2.0_f64*sll_pi))   .AND. &
        ( sqrt(normH1_0) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_1) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_2) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_3) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2)) then
     
     print *, 'PASSED'
  end if

contains

  ! obsolete apparently
  subroutine test_value_for_acceptable_error( value, max_error, boolean )
    sll_real64, intent(in) :: value
    sll_real64, intent(in) :: max_error
    logical, intent(inout) :: boolean

    if( value <= max_error ) then
       boolean = boolean .and. .true.
    else
       boolean = boolean .and. .false.
    end if
  end subroutine test_value_for_acceptable_error

end program unit_test

 
