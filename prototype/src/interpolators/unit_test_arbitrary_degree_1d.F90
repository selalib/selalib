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

  type(sll_arb_deg_1d_interpolator) :: ad1d
  sll_real64, dimension(:), allocatable    :: x,xprime
  sll_real64, dimension(:), allocatable    :: eta1_pos,eta1_prime
  sll_real64, dimension(:), allocatable    :: reference
  !sll_real64, dimension(:), allocatable    :: calculated
  !sll_real64, dimension(:), allocatable    :: difference
  !sll_int32 :: ierr
  sll_int32  :: i !, j
  sll_real64 :: eta1, h1
  sll_real64 :: acc, acc1, acc2,acc3,acc4,acc5,acc6,acc7
  sll_real64 :: node_val, ref, deriv1_val
  sll_real64 :: acc_der1,acc1_der1, acc2_der1,acc3_der1
  sll_real64 :: acc4_der1,acc5_der1,acc6_der1,acc7_der1
  sll_real64 :: normL2_0, normL2_1,normL2_2,normL2_3
  sll_real64 :: normL2_4,normL2_5,normL2_6,normL2_7
  sll_real64 :: normH1_0,normH1_1,normH1_2,normH1_3
  sll_real64 :: normH1_4,normH1_5,normH1_6,normH1_7
  
  print *,  'filling out discrete arrays for x1 '
  h1 = (X1MAX-X1MIN)/real(NPTS1-1,f64)
  print *, 'h1 = ', h1
  
  allocate(x(NPTS1))
  allocate(xprime(2))
  allocate(reference(NPTS1))
  !allocate(calculated(NPTS1))
  !allocate(difference(NPTS1))
  allocate(eta1_pos(NPTS1))
  allocate(eta1_prime(2))

  print *, '***********************************************************'
  print *, '              periodic case'
  print *, '***********************************************************'
  
  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = sin(2.0_f64*sll_pi*eta1)
     reference(i+1)     = sin(2.0_f64*sll_pi*eta1)
  end do

  
  
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
  normL2_0 = 0.0_f64
  normH1_0 = 0.0_f64
  do i=0,NPTS1-2
     eta1       = X1MIN + real(i,f64)*h1
     !print*, 'hi'
     node_val   = ad1d%interpolate_value(eta1)
     !print*, 'hi',node_val
     ref        = sin(2.0_f64*sll_pi*eta1)
    ! print*, 'hi',node_val, ref
     !calculated(i+1) = node_val
     !difference(i+1) = ref-node_val
     acc        = acc + abs(node_val-ref)
     normL2_0   = normL2_0  + (node_val-ref)**2*h1
     
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)
     ref = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
     ! print*, 'hi',deriv1_val, ref
     acc_der1 = acc_der1 + abs(deriv1_val-ref)
     normH1_0   = normH1_0  + (deriv1_val-ref)**2*h1
     !
  end do

!!$  call sll_gnuplot_1d(X1MIN, X1MAX, NPTS1, &
!!$       calculated, 'calculated_interp_arb_deg', 0, ierr)
!!$
!!$  call sll_gnuplot_1d(X1MIN, X1MAX, NPTS1, &
!!$       difference, 'difference_interp_arb_deg', 0, ierr)
!!$  
!!$  call sll_gnuplot_1d(X1MIN, X1MAX, NPTS1, &
!!$       ad1d%coeff_splines, 'coefficients_interp_arb_deg', 0, ierr)
  
  call delete(ad1d)
  
  print *, '***********************************************************'
  print *, '              dirichlet case'
  print *, '***********************************************************'
  

  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = sin(2.0_f64*sll_pi*eta1)
     reference(i+1)     = sin(2.0_f64*sll_pi*eta1)
  end do
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPL_DEG)
  
  call ad1d%compute_interpolants( &
       x(1:NPTS1))
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc1 = 0.0_f64
  acc1_der1 = 0.0_f64
  normL2_1 = 0.0_f64
  normH1_1 = 0.0_f64
  do i=0,NPTS1-2
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = sin(2.0_f64*sll_pi*eta1)
     !print*, 'hi',node_val, ref
     !calculated(i+1) = node_val
     !difference(i+1) = ref-node_val
     acc1        = acc1 + abs(node_val-ref)
     normL2_1   = normL2_1  + (node_val-ref)**2*h1
     !print*, 'dif', node_val-ref
     !print*, acc1
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
    ! print*, 'hi',deriv1_val, ref
     acc1_der1 = acc1_der1 + abs(deriv1_val-ref)
     normH1_1   = normH1_1  + (deriv1_val-ref)**2*h1
  end do
  
  call delete(ad1d)
  

  print *, '***********************************************************'
  print *, '              dirichlet case non homogene'
  print *, '***********************************************************'
  

  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = sin(2.0_f64*sll_pi*eta1) + 3.0_f64
     reference(i+1)     = sin(2.0_f64*sll_pi*eta1) + 3.0_f64
  end do
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPL_DEG)

  call set_values_at_boundary1d(ad1d,value_left=3.0_f64,value_right=3.0_f64)
  
  call ad1d%compute_interpolants( &
       x(1:NPTS1))
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc2 = 0.0_f64
  acc2_der1 = 0.0_f64
  normL2_2 = 0.0_f64
  normH1_2 = 0.0_f64
  do i=0,NPTS1-2
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = sin(2.0_f64*sll_pi*eta1)+ 3.0_f64
     !print*, 'hi',node_val, ref
     !calculated(i+1) = node_val
     !difference(i+1) = ref-node_val
     acc2        = acc2 + abs(node_val-ref)
     normL2_2    = normL2_2  + (node_val-ref)**2*h1
     !print*, 'dif', node_val-ref
     !print*, acc1
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref        = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
     acc2_der1  = acc2_der1 + abs(deriv1_val-ref)
     normH1_2   = normH1_2  + (deriv1_val-ref)**2*h1
  end do
  
  call delete(ad1d)


  print *, '***********************************************************'
  print *, '              Hermite-Hermite case '
  print *, '***********************************************************'
  

  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = sin(2.0_f64*sll_pi*eta1) + 3.0_f64
     reference(i+1)     = sin(2.0_f64*sll_pi*eta1) + 3.0_f64
  end do
  xprime(1) = cos(2.0_f64*sll_pi*eta1_pos(1))*2.0_f64*sll_pi
  xprime(2) = cos(2.0_f64*sll_pi*eta1_pos(NPTS1))*2.0_f64*sll_pi
  eta1_prime(1) = eta1_pos(1)
  eta1_prime(2) = eta1_pos(NPTS1)
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_HERMITE, &
       SLL_HERMITE, &
       SPL_DEG)

  call set_values_at_boundary1d(ad1d,value_left=3.0_f64,value_right=3.0_f64,&
       slope_left=xprime(1),slope_right=xprime(2))
  
  call ad1d%compute_interpolants(x(1:NPTS1))
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc3 = 0.0_f64
  acc3_der1 = 0.0_f64
  normL2_3 = 0.0_f64
  normH1_3 = 0.0_f64
  do i=0,NPTS1-2
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = sin(2.0_f64*sll_pi*eta1)+ 3.0_f64
     !print*, 'hi',node_val, ref
     !calculated(i+1) = node_val
     !difference(i+1) = ref-node_val
     acc3        = acc3 + abs(node_val-ref)
     normL2_3    = normL2_3  + (node_val-ref)**2*h1
     !print*, 'dif', node_val,ref
     !print*, acc1
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref        = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
     !print*, 'ho',deriv1_val, ref
     acc3_der1  = acc2_der1 + abs(deriv1_val-ref)
     normH1_3   = normH1_3  + (deriv1_val-ref)**2*h1
  end do
  
  call delete(ad1d)


  print *, '***********************************************************'
  print *, '              Neumann-Dirichlet case '
  print *, '***********************************************************'
  

  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = cos(2.0_f64*sll_pi*eta1)
     reference(i+1)     = cos(2.0_f64*sll_pi*eta1)
  end do
  xprime(1) = 0.0
  xprime(2) = 0.0
  eta1_prime(1) = eta1_pos(1)
  eta1_prime(2) = eta1_pos(NPTS1)
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_NEUMANN, &
       SLL_DIRICHLET, &
       SPL_DEG)

  call set_values_at_boundary1d(ad1d,value_left=1.0_f64,value_right=1.0_f64,&
       slope_left=xprime(1),slope_right=xprime(2))
  
  print*, 'ok'
  call ad1d%compute_interpolants(x(1:NPTS1))
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc4 = 0.0_f64
  acc4_der1 = 0.0_f64
  normL2_4 = 0.0_f64
  normH1_4 = 0.0_f64
  do i=0,NPTS1-2
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = cos(2.0_f64*sll_pi*eta1)
     !calculated(i+1) = node_val
     !difference(i+1) = ref-node_val
     acc4        = acc4 + abs(node_val-ref)
     normL2_4    = normL2_4  + (node_val-ref)**2*h1
     !print*, 'dif', node_val,ref
     !print*, acc1
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
     acc4_der1  = acc4_der1 + abs(deriv1_val-ref)
     normH1_4   = normH1_4  + (deriv1_val-ref)**2*h1
  end do
  
  call delete(ad1d)

  
  print *, '***********************************************************'
  print *, '              Neumann-Dirichlet case '
  print *, '***********************************************************'
  

  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = cos(2.0_f64*sll_pi*eta1)
     reference(i+1)     = cos(2.0_f64*sll_pi*eta1)
  end do
  xprime(1) = 0.0
  xprime(2) = 0.0
  eta1_prime(1) = eta1_pos(1)
  eta1_prime(2) = eta1_pos(NPTS1)
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_NEUMANN, &
       SLL_DIRICHLET, &
       SPL_DEG)

  call set_values_at_boundary1d(ad1d,value_left=1.0_f64,value_right=1.0_f64,&
       slope_left=xprime(1),slope_right=xprime(2))
  
  print*, 'ok'
  call ad1d%compute_interpolants(x(1:NPTS1))
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc4 = 0.0_f64
  acc4_der1 = 0.0_f64
  normL2_4 = 0.0_f64
  normH1_4 = 0.0_f64
  do i=0,NPTS1-1
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = cos(2.0_f64*sll_pi*eta1)
     acc4        = acc4 + abs(node_val-ref)
     normL2_4    = normL2_4  + (node_val-ref)**2*h1
    
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
     acc4_der1  = acc4_der1 + abs(deriv1_val-ref)
     normH1_4   = normH1_4  + (deriv1_val-ref)**2*h1
  end do
  
  call delete(ad1d)


  print *, '***********************************************************'
  print *, '              Neumann-Neumann case '
  print *, '***********************************************************'
  
  
  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = cos(2.0_f64*sll_pi*eta1)
     reference(i+1)     = cos(2.0_f64*sll_pi*eta1)
  end do
  xprime(1) = 0.0
  xprime(2) = 0.0
  eta1_prime(1) = eta1_pos(1)
  eta1_prime(2) = eta1_pos(NPTS1)
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_NEUMANN, &
       SLL_NEUMANN, &
       SPL_DEG)

  call set_values_at_boundary1d(ad1d,value_left=1.0_f64,value_right=1.0_f64,&
       slope_left=xprime(1),slope_right=xprime(2))
  
  print*, 'ok'
  call ad1d%compute_interpolants(x(1:NPTS1))
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc5 = 0.0_f64
  acc5_der1 = 0.0_f64
  normL2_5 = 0.0_f64
  normH1_5 = 0.0_f64
  do i=0,NPTS1-1
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = cos(2.0_f64*sll_pi*eta1)
     !print*, 'hi',node_val, ref
     !calculated(i+1) = node_val
     !difference(i+1) = ref-node_val
     acc5        = acc5 + abs(node_val-ref)
     normL2_5    = normL2_5  + (node_val-ref)**2*h1
     !print*, 'dif', node_val,ref
     !print*, acc1
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
     !print*, 'ho',deriv1_val, ref
     acc5_der1  = acc5_der1 + abs(deriv1_val-ref)
     normH1_5   = normH1_5  + (deriv1_val-ref)**2*h1
  end do
  
  call delete(ad1d)


  print *, '***********************************************************'
  print *, '              Hermite-Neumann case '
  print *, '***********************************************************'
  
  
  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = cos(2.0_f64*sll_pi*eta1)
     reference(i+1)     = cos(2.0_f64*sll_pi*eta1)
  end do
  xprime(1) = 0.0
  xprime(2) = 0.0
  eta1_prime(1) = eta1_pos(1)
  eta1_prime(2) = eta1_pos(NPTS1)
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_HERMITE, &
       SLL_NEUMANN, &
       SPL_DEG)

  call set_values_at_boundary1d(ad1d,value_left=1.0_f64,value_right=1.0_f64,&
       slope_left=xprime(1),slope_right=xprime(2))
  
  print*, 'ok'
  call ad1d%compute_interpolants(x(1:NPTS1))
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc6 = 0.0_f64
  acc6_der1 = 0.0_f64
  normL2_6 = 0.0_f64
  normH1_6 = 0.0_f64
  do i=0,NPTS1-1
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = cos(2.0_f64*sll_pi*eta1)
     !print*, 'hi',node_val, ref
     !calculated(i+1) = node_val
     !difference(i+1) = ref-node_val
     acc6        = acc6 + abs(node_val-ref)
     normL2_6    = normL2_6  + (node_val-ref)**2*h1
     !print*, 'dif', node_val,ref
     !print*, acc1
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
     !print*, 'ho',deriv1_val, ref
     acc6_der1  = acc6_der1 + abs(deriv1_val-ref)
     normH1_6   = normH1_6  + (deriv1_val-ref)**2*h1
  end do
  
  call delete(ad1d)
  
    print *, '***********************************************************'
  print *, '              Neumann-Hermite case '
  print *, '***********************************************************'
  
  
  do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = cos(2.0_f64*sll_pi*eta1)
     reference(i+1)     = cos(2.0_f64*sll_pi*eta1)
  end do
  xprime(1) = 0.0
  xprime(2) = 0.0
  eta1_prime(1) = eta1_pos(1)
  eta1_prime(2) = eta1_pos(NPTS1)
  
  call ad1d%initialize( &
       NPTS1, &
       X1MIN, &
       X1MAX, &
       SLL_NEUMANN, &
       SLL_HERMITE, &
       SPL_DEG)

  call set_values_at_boundary1d(ad1d,value_left=1.0_f64,value_right=1.0_f64,&
       slope_left=xprime(1),slope_right=xprime(2))
  
  print*, 'ok'
  call ad1d%compute_interpolants(x(1:NPTS1))
  
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc7 = 0.0_f64
  acc7_der1 = 0.0_f64
  normL2_7 = 0.0_f64
  normH1_7 = 0.0_f64
  do i=0,NPTS1-1
     eta1       = X1MIN + real(i,f64)*h1
     node_val   = ad1d%interpolate_value(eta1)
     ref        = cos(2.0_f64*sll_pi*eta1)
     !print*, 'hi',node_val, ref
     !calculated(i+1) = node_val
     !difference(i+1) = ref-node_val
     acc7        = acc7 + abs(node_val-ref)
     normL2_7    = normL2_7  + (node_val-ref)**2*h1
     !print*, 'dif', node_val,ref
     !print*, acc1
     deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
     ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
     !print*, 'ho',deriv1_val, ref
     acc7_der1  = acc7_der1 + abs(deriv1_val-ref)
     normH1_7   = normH1_7  + (deriv1_val-ref)**2*h1
  end do
  
  call delete(ad1d)


  print*, '--------------------------------------------'
  print*, ' Average error in nodes'
  print*, '--------------------------------------------'
  print *, 'Average error in nodes (periodic) = ', acc/(NPTS1)
  print *, 'Average error in nodes (dirichlet) = ', acc1/(NPTS1)
  print *, 'Average error in nodes (dirichlet non homogene) = ', acc2/(NPTS1)
  print *, 'Average error in nodes (Hermite-Hermite) = ', acc3/(NPTS1)
  print *, 'Average error in nodes (Neumann-Dirichlet) = ', acc4/(NPTS1)
  print *, 'Average error in nodes (Neumann-Neumann) = ', acc5/(NPTS1)
  print *, 'Average error in nodes (Hermite-Neumann) = ', acc6/(NPTS1)
  print *, 'Average error in nodes (Neumann-Hermite) = ', acc7/(NPTS1)
  print*, '--------------------------------------------'
  print*, ' Average error in nodes first derivative eta1'
  print*, '--------------------------------------------'
 
  print *,'Average error in nodes first derivative eta1(periodic)=',&
       acc_der1/(NPTS1)
  print *,'Average error in nodes first derivative eta1(dirichlet)=',&
       acc1_der1/(NPTS1)
  print *,'Average error in nodes first derivative eta1(dirichlet non homogene)=',&
       acc2_der1/(NPTS1)
  print *,'Average error in nodes first derivative eta1(Hermite-Hermite)=',&
       acc3_der1/(NPTS1)
  print *,'Average error in nodes first derivative eta1(Neumann-Dirichlet)=',&
       acc4_der1/(NPTS1)
  print *,'Average error in nodes first derivative eta1(Neumann-Neumann)=',&
       acc5_der1/(NPTS1)
  print *,'Average error in nodes first derivative eta1(Hermite-Neumann)=',&
       acc6_der1/(NPTS1)
  print *,'Average error in nodes first derivative eta1(Neumann-Hermite)=',&
       acc7_der1/(NPTS1)
 
  print*, '--------------------------------------------'
  print*, ' Norm L2 error '
  print*, '--------------------------------------------'
  print*,  'norm L2 error (periodic) =', sqrt(normL2_0), h1**(SPL_DEG)
  print*,  'norm L2 error (dirichlet) =', sqrt(normL2_1),h1**(SPL_DEG)
  print*,  'norm L2 error (dirichlet non homogene) =', sqrt(normL2_2),h1**(SPL_DEG)
  print*,  'norm L2 error (Hermite-Hermite) =', sqrt(normL2_3),h1**(SPL_DEG)
  print*,  'norm L2 error (Neumann-dirichlet) =', sqrt(normL2_4),h1**(SPL_DEG)
  print*,  'norm L2 error (Neumann-Neumann) =', sqrt(normL2_5),h1**(SPL_DEG)
  print*,  'norm L2 error (Hermite-Neumann) =', sqrt(normL2_6),h1**(SPL_DEG)
  print*,  'norm L2 error (Neumann-Hermite) =', sqrt(normL2_7),h1**(SPL_DEG)

  print*, '--------------------------------------------'
  print*, ' Norm H1 error '
  print*, '--------------------------------------------'
  print*,  'norm H1 error (periodic) =',  sqrt(normH1_0),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
  print*,  'norm H1 error (dirichlet) =', sqrt(normH1_1),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
  print*,  'norm H1 error (dirichlet non homogene) =', sqrt(normH1_2),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
  print*,  'norm H1 error (Hermite-Hermite) =',sqrt(normH1_3),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
  print*,  'norm H1 error (Neumann-Dirichlet) =',sqrt(normH1_4),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
  print*,  'norm H1 error (Neumann-Neumann) =',sqrt(normH1_5),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
  print*,  'norm H1 error (Hermite-Neumann) =',sqrt(normH1_6),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
  print*,  'norm H1 error (Neumann-Hermite) =',sqrt(normH1_7),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2

  if (  ( sqrt(normL2_0) <= h1**(SPL_DEG+1)) .AND. &
        ( sqrt(normL2_1) <= h1**(SPL_DEG+1)) .AND. &
        ( sqrt(normL2_2) <= h1**(SPL_DEG+1)) .AND. & 
        ( sqrt(normL2_3) <= h1**(SPL_DEG+1)) .AND. &
        ( sqrt(normL2_4) <= h1**(SPL_DEG+1)) .AND. &
        ( sqrt(normL2_5) <= h1**(SPL_DEG+1)) .AND. &
        ( sqrt(normL2_6) <= h1**(SPL_DEG+1)) .AND. &
        ( sqrt(normL2_7) <= h1**(SPL_DEG+1)) .AND. &
        ( sqrt(normH1_0) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_1) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_2) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_3) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_4) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_5) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_6) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) .AND. &
        ( sqrt(normH1_7) <= h1**(SPL_DEG-3)*(2.0_f64*sll_pi)**2) ) then
     
       
     print *, 'PASSED'
  end if

  deallocate(x)
  deallocate(xprime)
  deallocate(reference)
  !deallocate(calculated)
  !deallocate(difference)
  deallocate(eta1_pos)
  deallocate(eta1_prime)

end program unit_test

 
