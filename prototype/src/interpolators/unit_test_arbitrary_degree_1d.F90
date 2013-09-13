program unit_test
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_interpolators.h"
  use sll_arbitrary_degree_spline_interpolator_1d_module
  use sll_gnuplot
  implicit none

#define NPTS1 31
#define NPTS2 31 
#define SPL_DEG 3
#define X1MIN 0.0_f64
#define X1MAX 1.0_f64

  type(arb_deg_1d_interpolator) :: ad1d
  sll_real64, dimension(:), allocatable    :: x
  sll_real64, dimension(:), allocatable      :: eta1_pos
  sll_real64, dimension(:), allocatable    :: reference
  sll_real64, dimension(:), allocatable    :: calculated
  sll_real64, dimension(:), allocatable    :: difference
  sll_int32 :: ierr
  sll_int32  :: i, j
  sll_real64 :: eta1, h1
  sll_real64 :: acc, acc1, acc2, acc3, node_val, ref, deriv1_val
  sll_real64 :: acc_der1, acc1_der1, acc2_der1, acc3_der1

  
  print *,  'filling out discrete arrays for x1 '
  h1 = (X1MAX-X1MIN)/real(NPTS1-1,f64)
  print *, 'h1 = ', h1
  
  allocate(x(NPTS1))
  allocate(reference(NPTS1))
  allocate(calculated(NPTS1))
  allocate(difference(NPTS1))
  allocate(eta1_pos(NPTS1))

  print *, '***********************************************************'
  print *, '              periodic case'
  print *, '***********************************************************'
  
  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = sin(2.0_f64*sll_pi*eta1)
     reference(i+1)     = sin(2.0_f64*sll_pi*eta1)
  end do
!!$  
!!$  call sll_gnuplot_field_1d( &
!!$       X1MIN, &
!!$       NPTS1, &
!!$       reference, &
!!$       'reference_interp_arb_deg', &
!!$       0, &
!!$       ierr)
  
  
  ! Test the 1D transformation:
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SPL_DEG)
  
  call ad1d%compute_interpolants( &
       x(1:NPTS1-1))!,&
       !eta1_pos(1:NPTS1-1),&
       !NPTS1-1)
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  acc  = 0.0_f64
  acc_der1  = 0.0_f64
  
  do i=0,NPTS1-2
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = sin(2.0_f64*sll_pi*eta1)
     calculated(i+1) = node_val
     difference(i+1) = ref-node_val
     acc        = acc + abs(node_val-ref)
     
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)
     ref = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
     acc_der1 = acc_der1 + abs(deriv1_val-ref)
     !
  end do

!!$  call sll_gnuplot_field_1d(X1MIN, X1MAX, NPTS1, &
!!$       calculated, 'calculated_interp_arb_deg', 0, ierr)
!!$
!!$  call sll_gnuplot_field_1d(X1MIN, X1MAX, NPTS1, &
!!$       difference, 'difference_interp_arb_deg', 0, ierr)
!!$  
!!$  call sll_gnuplot_field_1d(X1MIN, X1MAX, NPTS1, &
!!$       ad1d%coeff_splines, 'coefficients_interp_arb_deg', 0, ierr)
  
  
  print *, '***********************************************************'
  print *, '              dirichlet case'
  print *, '***********************************************************'
  
  call delete(ad1d)
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPL_DEG)
  
  call ad1d%compute_interpolants( &
       x(1:NPTS1-1))!,&
       !eta1_pos(1:NPTS1-1),&
       !NPTS1-1)
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc1 = 0.0_f64
  acc1_der1 = 0.0_f64
  
  do i=0,NPTS1-2
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = sin(2.0_f64*sll_pi*eta1)
     calculated(i+1) = node_val
     difference(i+1) = ref-node_val
     acc1        = acc1 + abs(node_val-ref)
     
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
     acc1_der1 = acc1_der1 + abs(deriv1_val-ref)
  end do
  
  call delete(ad1d)
  
  print*, '--------------------------------------------'
  print*, ' Average error in nodes'
  print*, '--------------------------------------------'
  print *, 'Average error in nodes (dirichlet) = ', acc1/(NPTS1)
!
  print *, 'Average error in nodes (periodic) = ', acc/(NPTS1)

  print*, '--------------------------------------------'
  print*, ' Average error in nodes first derivative eta1'
  print*, '--------------------------------------------'
 
  print *,'Average error in nodes first derivative eta1(dirichlet)=',&
       acc1_der1/(NPTS1)
 
  print *,'Average error in nodes first derivative eta1(periodic)=',&
       acc_der1/(NPTS1)

 
!!$
!!$  if( (acc/(NPTS1*NPTS2)  .lt. 2.0e-16) .and. &
!!$      (acc1/(NPTS1*NPTS2) .lt. 2.0e-16) .and. &
!!$     print *, 'PASSED'
!!$  else
!!$     print *, 'FAILED'
!!$  end if
!!$  print *, 'Average error, x1 deriv eta1 = ', acc1/(NPTS1)
!!$  print *, 'Average error, x1 deriv eta1 = ', acc/(NPTS1)

end program unit_test

 
