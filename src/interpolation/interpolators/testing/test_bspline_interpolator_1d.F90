program unit_test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_bspline_interpolator_1d, only: &
    sll_s_set_values_at_boundary1d, &
    sll_t_bspline_interpolator_1d, &
    sll_o_delete

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_hermite, &
    sll_p_neumann, &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NPTS1 31
#define NPTS2 31 
#define SPL_DEG 3
#define X1MIN 0.0_f64
#define X1MAX 1.0_f64

type(sll_t_bspline_interpolator_1d) :: ad1d

sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: xprime
sll_real64, dimension(:), allocatable :: eta1_pos
sll_real64, dimension(:), allocatable :: eta1_prime
sll_real64, dimension(:), allocatable :: reference

sll_int32  :: i 
sll_real64 :: eta1, h1
sll_real64 :: acc, acc1, acc2,acc3,acc4,acc5,acc6,acc7
sll_real64 :: node_val, ref, deriv1_val
sll_real64 :: acc_der1,acc1_der1, acc2_der1,acc3_der1
sll_real64 :: acc4_der1,acc5_der1,acc6_der1,acc7_der1
sll_real64 :: normL2_0, normL2_1,normL2_2,normL2_3
sll_real64 :: normL2_4,normL2_5,normL2_6,normL2_7
sll_real64 :: normH1_0,normH1_1,normH1_2,normH1_3
sll_real64 :: normH1_4,normH1_5,normH1_6,normH1_7
sll_real64, parameter :: dpi   = 2*sll_p_pi
sll_real64, parameter :: dpisq = dpi*dpi
  
h1 = (X1MAX-X1MIN)/real(NPTS1-1,f64)
  
allocate(x(NPTS1))
allocate(xprime(2))
allocate(reference(NPTS1))
allocate(eta1_pos(NPTS1))
allocate(eta1_prime(2))

do i=0,NPTS1-1
  eta1           = X1MIN + real(i,f64)*h1
  eta1_pos(i+1)  = eta1
  x(i+1)         = sin(2.0_f64*sll_p_pi*eta1)
  reference(i+1) = sin(2.0_f64*sll_p_pi*eta1)
end do
  
call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_periodic,sll_p_periodic,SPL_DEG)
call ad1d%compute_interpolants(x)
  
acc  = 0.0_f64
acc_der1  = 0.0_f64
normL2_0 = 0.0_f64
normH1_0 = 0.0_f64
do i=0,NPTS1-2
  eta1       = X1MIN + real(i,f64)*h1
  node_val   = ad1d%interpolate_from_interpolant_value(eta1)
  ref        = sin(2.0_f64*sll_p_pi*eta1)
  acc        = acc + abs(node_val-ref)
  normL2_0   = normL2_0  + (node_val-ref)**2*h1
  deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)
  ref        = 2.0_f64*sll_p_pi*cos(2.0_f64*sll_p_pi*eta1)
  acc_der1   = acc_der1 + abs(deriv1_val-ref)
  normH1_0   = normH1_0  + (deriv1_val-ref)**2*h1
end do

call sll_o_delete(ad1d)
  
do i=0,NPTS1-1
  eta1           = X1MIN + real(i,f64)*h1
  eta1_pos(i+1)  = eta1
  x(i+1)         = sin(2.0_f64*sll_p_pi*eta1)
  reference(i+1) = sin(2.0_f64*sll_p_pi*eta1)
end do
  
call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_dirichlet,sll_p_dirichlet,SPL_DEG)
call ad1d%compute_interpolants(x)
  
acc1 = 0.0_f64
acc1_der1 = 0.0_f64
normL2_1 = 0.0_f64
normH1_1 = 0.0_f64
do i=0,NPTS1-2
   eta1       = X1MIN + real(i,f64)*h1
   node_val   = ad1d%interpolate_from_interpolant_value(eta1)
   ref        = sin(2.0_f64*sll_p_pi*eta1)
   acc1       = acc1 + abs(node_val-ref)
   normL2_1   = normL2_1  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)   
   ref        = 2.0_f64*sll_p_pi*cos(2.0_f64*sll_p_pi*eta1)
   acc1_der1  = acc1_der1 + abs(deriv1_val-ref)
   normH1_1   = normH1_1  + (deriv1_val-ref)**2*h1
end do

call sll_o_delete(ad1d)
  
do i=0,NPTS1-1
  eta1           = X1MIN + real(i,f64)*h1
  eta1_pos(i+1)  = eta1
  x(i+1)         = sin(2.0_f64*sll_p_pi*eta1) + 3.0_f64
  reference(i+1) = sin(2.0_f64*sll_p_pi*eta1) + 3.0_f64
end do
  
call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_dirichlet,sll_p_dirichlet,SPL_DEG)

call sll_s_set_values_at_boundary1d(ad1d,value_left=3.0_f64,value_right=3.0_f64)
  
call ad1d%compute_interpolants(x)
  
acc2 = 0.0_f64
acc2_der1 = 0.0_f64
normL2_2 = 0.0_f64
normH1_2 = 0.0_f64
do i=0,NPTS1-2
  eta1       = X1MIN + real(i,f64)*h1
  node_val   = ad1d%interpolate_from_interpolant_value(eta1)
  ref        = sin(2.0_f64*sll_p_pi*eta1)+ 3.0_f64
  acc2       = acc2 + abs(node_val-ref)
  normL2_2   = normL2_2  + (node_val-ref)**2*h1
  deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)   
  ref        = 2.0_f64*sll_p_pi*cos(2.0_f64*sll_p_pi*eta1)
  acc2_der1  = acc2_der1 + abs(deriv1_val-ref)
  normH1_2   = normH1_2  + (deriv1_val-ref)**2*h1
end do
  
call sll_o_delete(ad1d)

do i=0,NPTS1-1
     eta1               = X1MIN + real(i,f64)*h1
     eta1_pos(i+1)      = eta1
     x(i+1)             = sin(2.0_f64*sll_p_pi*eta1) + 3.0_f64
     reference(i+1)     = sin(2.0_f64*sll_p_pi*eta1) + 3.0_f64
end do
xprime(1)     = cos(2.0_f64*sll_p_pi*eta1_pos(1))*2.0_f64*sll_p_pi
xprime(2)     = cos(2.0_f64*sll_p_pi*eta1_pos(NPTS1))*2.0_f64*sll_p_pi
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)
  
call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_hermite,sll_p_hermite,SPL_DEG)

call sll_s_set_values_at_boundary1d(ad1d,                   &
                              value_left=3.0_f64,     &
                              value_right=3.0_f64,    &
                              slope_left=xprime(1),   &
                              slope_right=xprime(2))
 
call ad1d%compute_interpolants(x)
  
acc3 = 0.0_f64
acc3_der1 = 0.0_f64
normL2_3 = 0.0_f64
normH1_3 = 0.0_f64
do i=0,NPTS1-2
  eta1       = X1MIN + real(i,f64)*h1
  node_val   = ad1d%interpolate_from_interpolant_value(eta1)
  ref        = sin(2.0_f64*sll_p_pi*eta1)+ 3.0_f64
  acc3       = acc3 + abs(node_val-ref)
  normL2_3   = normL2_3  + (node_val-ref)**2*h1
  deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)   
  ref        = 2.0_f64*sll_p_pi*cos(2.0_f64*sll_p_pi*eta1)
  acc3_der1  = acc2_der1 + abs(deriv1_val-ref)
  normH1_3   = normH1_3  + (deriv1_val-ref)**2*h1
end do
  
call sll_o_delete(ad1d)

do i=0,NPTS1-1
  eta1               = X1MIN + real(i,f64)*h1
  eta1_pos(i+1)      = eta1
  x(i+1)             = cos(2.0_f64*sll_p_pi*eta1)
  reference(i+1)     = cos(2.0_f64*sll_p_pi*eta1)
end do

xprime(1) = 0.0_f64
xprime(2) = 0.0_f64
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_neumann,sll_p_dirichlet,SPL_DEG)

call sll_s_set_values_at_boundary1d(ad1d,                 &
                              value_left=1.0_f64,   &
                              value_right=1.0_f64,  &
                              slope_left=xprime(1), &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc4 = 0.0_f64
acc4_der1 = 0.0_f64
normL2_4 = 0.0_f64
normH1_4 = 0.0_f64
do i=0,NPTS1-2
   eta1       = X1MIN + real(i,f64)*h1
   node_val   = ad1d%interpolate_from_interpolant_value(eta1)
   ref        = cos(2.0_f64*sll_p_pi*eta1)
   acc4       = acc4 + abs(node_val-ref)
   normL2_4   = normL2_4  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_p_pi*sin(2.0_f64*sll_p_pi*eta1)
   acc4_der1  = acc4_der1 + abs(deriv1_val-ref)
   normH1_4   = normH1_4  + (deriv1_val-ref)**2*h1
end do

call sll_o_delete(ad1d)

do i=0,NPTS1-1
   eta1           = X1MIN + real(i,f64)*h1
   eta1_pos(i+1)  = eta1
   x(i+1)         = cos(2.0_f64*sll_p_pi*eta1)
   reference(i+1) = cos(2.0_f64*sll_p_pi*eta1)
end do
xprime(1) = 0.0_f64
xprime(2) = 0.0_f64
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_neumann,sll_p_dirichlet,SPL_DEG)

call sll_s_set_values_at_boundary1d(ad1d,                 &
                              value_left=1.0_f64,   &
                              value_right=1.0_f64,  &
                              slope_left=xprime(1), &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc4 = 0.0_f64
acc4_der1 = 0.0_f64
normL2_4 = 0.0_f64
normH1_4 = 0.0_f64
do i=0,NPTS1-1
   eta1       = X1MIN + real(i,f64)*h1
   node_val   = ad1d%interpolate_from_interpolant_value(eta1)
   ref        = cos(2.0_f64*sll_p_pi*eta1)
   acc4       = acc4 + abs(node_val-ref)
   normL2_4   = normL2_4  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_p_pi*sin(2.0_f64*sll_p_pi*eta1)
   acc4_der1  = acc4_der1 + abs(deriv1_val-ref)
   normH1_4   = normH1_4  + (deriv1_val-ref)**2*h1
end do

call sll_o_delete(ad1d)

do i=0,NPTS1-1
   eta1           = X1MIN + real(i,f64)*h1
   eta1_pos(i+1)  = eta1
   x(i+1)         = cos(2.0_f64*sll_p_pi*eta1)
   reference(i+1) = cos(2.0_f64*sll_p_pi*eta1)
end do
xprime(1) = 0.0_f64
xprime(2) = 0.0_f64
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_neumann,sll_p_neumann,SPL_DEG)

call sll_s_set_values_at_boundary1d(ad1d,                 &
                              value_left=1.0_f64,   &
                              value_right=1.0_f64,  &
                              slope_left=xprime(1), &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc5 = 0.0_f64
acc5_der1 = 0.0_f64
normL2_5 = 0.0_f64
normH1_5 = 0.0_f64
do i=0,NPTS1-1
   eta1       = X1MIN + real(i,f64)*h1
   node_val   = ad1d%interpolate_from_interpolant_value(eta1)
   ref        = cos(2.0_f64*sll_p_pi*eta1)
   acc5       = acc5 + abs(node_val-ref)
   normL2_5   = normL2_5  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_p_pi*sin(2.0_f64*sll_p_pi*eta1)
   acc5_der1  = acc5_der1 + abs(deriv1_val-ref)
   normH1_5   = normH1_5  + (deriv1_val-ref)**2*h1
end do

call sll_o_delete(ad1d)

do i=0,NPTS1-1
   eta1           = X1MIN + real(i,f64)*h1
   eta1_pos(i+1)  = eta1
   x(i+1)         = cos(2.0_f64*sll_p_pi*eta1)
   reference(i+1) = cos(2.0_f64*sll_p_pi*eta1)
end do
xprime(1) = 0.0_f64
xprime(2) = 0.0_f64
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_hermite,sll_p_neumann,SPL_DEG)

call sll_s_set_values_at_boundary1d(ad1d,                  &
                              value_left=1.0_f64,    &
                              value_right=1.0_f64,   &
                              slope_left=xprime(1),  &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc6 = 0.0_f64
acc6_der1 = 0.0_f64
normL2_6 = 0.0_f64
normH1_6 = 0.0_f64
do i=0,NPTS1-1
   eta1       = X1MIN + real(i,f64)*h1
   node_val   = ad1d%interpolate_from_interpolant_value(eta1)
   ref        = cos(2.0_f64*sll_p_pi*eta1)
   acc6       = acc6 + abs(node_val-ref)
   normL2_6   = normL2_6  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_p_pi*sin(2.0_f64*sll_p_pi*eta1)
   acc6_der1  = acc6_der1 + abs(deriv1_val-ref)
   normH1_6   = normH1_6  + (deriv1_val-ref)**2*h1
end do

call sll_o_delete(ad1d)

do i=0,NPTS1-1
   eta1           = X1MIN + real(i,f64)*h1
   eta1_pos(i+1)  = eta1
   x(i+1)         = cos(2.0_f64*sll_p_pi*eta1)
   reference(i+1) = cos(2.0_f64*sll_p_pi*eta1)
end do
xprime(1) = 0.0_f64
xprime(2) = 0.0_f64
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,sll_p_neumann,sll_p_hermite,SPL_DEG)

call sll_s_set_values_at_boundary1d(ad1d,                 &
                              value_left=1.0_f64,   &
                              value_right=1.0_f64,  &
                              slope_left=xprime(1), &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc7 = 0.0_f64
acc7_der1 = 0.0_f64
normL2_7 = 0.0_f64
normH1_7 = 0.0_f64
do i=0,NPTS1-1
   eta1       = X1MIN + real(i,f64)*h1
   node_val   = ad1d%interpolate_from_interpolant_value(eta1)
   ref        = cos(2.0_f64*sll_p_pi*eta1)
   acc7        = acc7 + abs(node_val-ref)
   normL2_7    = normL2_7  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_from_interpolant_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_p_pi*sin(2.0_f64*sll_p_pi*eta1)
   acc7_der1  = acc7_der1 + abs(deriv1_val-ref)
   normH1_7   = normH1_7  + (deriv1_val-ref)**2*h1
end do

call sll_o_delete(ad1d)


write(*,100)'--------------------------------------------'
write(*,100)' Average error in nodes'
write(*,100)'--------------------------------------------'
write(*,100)'periodic               = ', acc/ real(NPTS1,f64)
write(*,100)'dirichlet              = ', acc1/real(NPTS1,f64)
write(*,100)'dirichlet non homogene = ', acc2/real(NPTS1,f64)
write(*,100)'Hermite-Hermite        = ', acc3/real(NPTS1,f64)
write(*,100)'Neumann-Dirichlet      = ', acc4/real(NPTS1,f64)
write(*,100)'Neumann-Neumann        = ', acc5/real(NPTS1,f64)
write(*,100)'Hermite-Neumann        = ', acc6/real(NPTS1,f64)
write(*,100)'Neumann-Hermite        = ', acc7/real(NPTS1,f64)
write(*,100)' '
write(*,100)'--------------------------------------------'
write(*,100)' Average error in nodes first derivative'
write(*,100)'--------------------------------------------'
write(*,100)'periodic               = ', acc_der1 /real(NPTS1,f64)
write(*,100)'dirichlet              = ', acc1_der1/real(NPTS1,f64)
write(*,100)'dirichlet non homogene = ', acc2_der1/real(NPTS1,f64)
write(*,100)'Hermite-Hermite        = ', acc3_der1/real(NPTS1,f64)
write(*,100)'Neumann-Dirichlet      = ', acc4_der1/real(NPTS1,f64)
write(*,100)'Neumann-Neumann        = ', acc5_der1/real(NPTS1,f64)
write(*,100)'Hermite-Neumann        = ', acc6_der1/real(NPTS1,f64)
write(*,100)'Neumann-Hermite        = ', acc7_der1/real(NPTS1,f64)
write(*,100) ' '
write(*,100)'--------------------------------------------'
write(*,100)' Norm L2 error '
write(*,100)'--------------------------------------------'
write(*,100)'periodic               = ', sqrt(normL2_0), h1**(SPL_DEG)
write(*,100)'dirichlet              = ', sqrt(normL2_1), h1**(SPL_DEG)
write(*,100)'dirichlet non homogene = ', sqrt(normL2_2), h1**(SPL_DEG)
write(*,100)'Hermite-Hermite        = ', sqrt(normL2_3), h1**(SPL_DEG)
write(*,100)'Neumann-dirichlet      = ', sqrt(normL2_4), h1**(SPL_DEG)
write(*,100)'Neumann-Neumann        = ', sqrt(normL2_5), h1**(SPL_DEG)
write(*,100)'Hermite-Neumann        = ', sqrt(normL2_6), h1**(SPL_DEG)
write(*,100)'Neumann-Hermite        = ', sqrt(normL2_7), h1**(SPL_DEG)
write(*,100) ' '
write(*,100)'--------------------------------------------'
write(*,100)' Norm H1 error '
write(*,100)'--------------------------------------------'
write(*,100)'periodic               = ', sqrt(normH1_0), h1**(SPL_DEG-2)*dpisq
write(*,100)'dirichlet              = ', sqrt(normH1_1), h1**(SPL_DEG-2)*dpisq
write(*,100)'dirichlet non homogene = ', sqrt(normH1_2), h1**(SPL_DEG-2)*dpisq
write(*,100)'Hermite-Hermite        = ', sqrt(normH1_3), h1**(SPL_DEG-2)*dpisq
write(*,100)'Neumann-Dirichlet      = ', sqrt(normH1_4), h1**(SPL_DEG-2)*dpisq
write(*,100)'Neumann-Neumann        = ', sqrt(normH1_5), h1**(SPL_DEG-2)*dpisq
write(*,100)'Hermite-Neumann        = ', sqrt(normH1_6), h1**(SPL_DEG-2)*dpisq
write(*,100)'Neumann-Hermite        = ', sqrt(normH1_7), h1**(SPL_DEG-2)*dpisq

100 format(a,2f25.20)

if(( sqrt(normL2_0) <= h1**(SPL_DEG+1)) .AND. &
   ( sqrt(normL2_1) <= h1**(SPL_DEG+1)) .AND. &
   ( sqrt(normL2_2) <= h1**(SPL_DEG+1)) .AND. & 
   ( sqrt(normL2_3) <= h1**(SPL_DEG+1)) .AND. &
   ( sqrt(normL2_4) <= h1**(SPL_DEG+1)) .AND. &
   ( sqrt(normL2_5) <= h1**(SPL_DEG+1)) .AND. &
   ( sqrt(normL2_6) <= h1**(SPL_DEG+1)) .AND. &
   ( sqrt(normL2_7) <= h1**(SPL_DEG+1)) .AND. &
   ( sqrt(normH1_0) <= h1**(SPL_DEG-3)*(2.0_f64*sll_p_pi)**2) .AND. &
   ( sqrt(normH1_1) <= h1**(SPL_DEG-3)*(2.0_f64*sll_p_pi)**2) .AND. &
   ( sqrt(normH1_2) <= h1**(SPL_DEG-3)*(2.0_f64*sll_p_pi)**2) .AND. &
   ( sqrt(normH1_3) <= h1**(SPL_DEG-3)*(2.0_f64*sll_p_pi)**2) .AND. &
   ( sqrt(normH1_4) <= h1**(SPL_DEG-3)*(2.0_f64*sll_p_pi)**2) .AND. &
   ( sqrt(normH1_5) <= h1**(SPL_DEG-3)*(2.0_f64*sll_p_pi)**2) .AND. &
   ( sqrt(normH1_6) <= h1**(SPL_DEG-3)*(2.0_f64*sll_p_pi)**2) .AND. &
   ( sqrt(normH1_7) <= h1**(SPL_DEG-3)*(2.0_f64*sll_p_pi)**2) ) then

  print *, 'PASSED'

end if

deallocate(x)
deallocate(xprime)
deallocate(reference)
deallocate(eta1_pos)
deallocate(eta1_prime)

end program unit_test
