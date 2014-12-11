program unit_test
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_interpolators.h"
implicit none

#define NPTS1 65
#define NPTS2 65 
#define SPL_DEG 3
#define X1MIN (-5.0_f64)
#define X1MAX (+5.0_f64)

type(sll_arbitrary_degree_spline_interpolator_1d) :: ad1d

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
  
h1 = (X1MAX-X1MIN)/real(NPTS1-1,f64)
  
allocate(x(NPTS1))
allocate(xprime(2))
allocate(reference(NPTS1))
allocate(eta1_pos(NPTS1))
allocate(eta1_prime(2))

!print *, '***********************************************************'
!print *, '              periodic case'
!print *, '***********************************************************'
  
do i=1,NPTS1
  eta1         = X1MIN + (i-1)*h1
  eta1_pos(i)  = eta1
  x(i)         = sin(2.0_f64*sll_pi*eta1)
  reference(i) = sin(2.0_f64*sll_pi*eta1)
end do
  
call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_PERIODIC,SLL_PERIODIC,SPL_DEG)
call ad1d%compute_interpolants(x)
  
acc  = 0.0_f64
acc_der1  = 0.0_f64
normL2_0 = 0.0_f64
normH1_0 = 0.0_f64
do i=1,NPTS1
  eta1       = eta1_pos(i)
  node_val   = ad1d%interpolate_value(eta1)
  ref        = reference(i)
  acc        = acc + abs(node_val-ref)
  normL2_0   = normL2_0  + (node_val-ref)**2*h1
  deriv1_val = ad1d%interpolate_derivative_eta1(eta1)
  ref        = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
  acc_der1   = acc_der1 + abs(deriv1_val-ref)
  normH1_0   = normH1_0  + (deriv1_val-ref)**2*h1
end do

call sll_delete(ad1d)
  
!print *, '***********************************************************'
!print *, '              dirichlet case'
!print *, '***********************************************************'
  
call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_DIRICHLET,SLL_DIRICHLET,SPL_DEG)
call ad1d%compute_interpolants(x)
  
acc1 = 0.0_f64
acc1_der1 = 0.0_f64
normL2_1 = 0.0_f64
normH1_1 = 0.0_f64
do i=1,NPTS1
   eta1       = eta1_pos(i)
   node_val   = ad1d%interpolate_value(eta1)
   ref        = reference(i)
   acc1       = acc1 + abs(node_val-ref)
   normL2_1   = normL2_1  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
   ref        = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
   acc1_der1  = acc1_der1 + abs(deriv1_val-ref)
   normH1_1   = normH1_1  + (deriv1_val-ref)**2*h1
end do

call sll_delete(ad1d)
  
!print *, '***********************************************************'
!print *, '              dirichlet case non homogene'
!print *, '***********************************************************'

do i=1,NPTS1
  eta1         = eta1_pos(i)
  x(i)         = sin(2.0_f64*sll_pi*eta1) + 3.0_f64
  reference(i) = sin(2.0_f64*sll_pi*eta1) + 3.0_f64
end do
  
call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_DIRICHLET,SLL_DIRICHLET,SPL_DEG)

call set_values_at_boundary1d(ad1d,value_left=3.0_f64,value_right=3.0_f64)
  
call ad1d%compute_interpolants(x)
  
acc2 = 0.0_f64
acc2_der1 = 0.0_f64
normL2_2 = 0.0_f64
normH1_2 = 0.0_f64
do i=1,NPTS1
  eta1       = eta1_pos(i)
  node_val   = ad1d%interpolate_value(eta1)
  ref        = sin(2.0_f64*sll_pi*eta1)+ 3.0_f64
  acc2       = acc2 + abs(node_val-ref)
  normL2_2   = normL2_2  + (node_val-ref)**2*h1
  deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
  ref        = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
  acc2_der1  = acc2_der1 + abs(deriv1_val-ref)
  normH1_2   = normH1_2  + (deriv1_val-ref)**2*h1
end do
  
call sll_delete(ad1d)

!print *, '***********************************************************'
!print *, '              Hermite-Hermite case '
!print *, '***********************************************************'

do i=1,NPTS1
  eta1         = eta1_pos(i)
  x(i)         = sin(2.0_f64*sll_pi*eta1) + 3.0_f64
  reference(i) = sin(2.0_f64*sll_pi*eta1) + 3.0_f64
end do
xprime(1)     = cos(2.0_f64*sll_pi*eta1_pos(1))*2.0_f64*sll_pi
xprime(2)     = cos(2.0_f64*sll_pi*eta1_pos(NPTS1))*2.0_f64*sll_pi
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)
  
call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_HERMITE,SLL_HERMITE,SPL_DEG)

call set_values_at_boundary1d(ad1d,                   &
                              value_left=3.0_f64,     &
                              value_right=3.0_f64,    &
                              slope_left=xprime(1),   &
                              slope_right=xprime(2))
 
call ad1d%compute_interpolants(x)
  
acc3 = 0.0_f64
acc3_der1 = 0.0_f64
normL2_3 = 0.0_f64
normH1_3 = 0.0_f64
do i=1,NPTS1
  eta1       = eta1_pos(i)
  node_val   = ad1d%interpolate_value(eta1)
  ref        = sin(2.0_f64*sll_pi*eta1)+ 3.0_f64
  acc3       = acc3 + abs(node_val-ref)
  normL2_3   = normL2_3  + (node_val-ref)**2*h1
  deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
  ref        = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*eta1)
  acc3_der1  = acc2_der1 + abs(deriv1_val-ref)
  normH1_3   = normH1_3  + (deriv1_val-ref)**2*h1
end do
  
call sll_delete(ad1d)

!print *, '***********************************************************'
!print *, '              Neumann-Dirichlet case '
!print *, '***********************************************************'

do i=1,NPTS1
  eta1         = eta1_pos(i)
  x(i)         = cos(2.0_f64*sll_pi*eta1)
  reference(i) = cos(2.0_f64*sll_pi*eta1)
end do

xprime(1) = 0.0
xprime(2) = 0.0
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_NEUMANN,SLL_DIRICHLET,SPL_DEG)

call set_values_at_boundary1d(ad1d,                 &
                              value_left=1.0_f64,   &
                              value_right=1.0_f64,  &
                              slope_left=xprime(1), &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc4 = 0.0_f64
acc4_der1 = 0.0_f64
normL2_4 = 0.0_f64
normH1_4 = 0.0_f64
do i=1,NPTS1
   eta1       = eta1_pos(i)
   node_val   = ad1d%interpolate_value(eta1)
   ref        = cos(2.0_f64*sll_pi*eta1)
   acc4       = acc4 + abs(node_val-ref)
   normL2_4   = normL2_4  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
   acc4_der1  = acc4_der1 + abs(deriv1_val-ref)
   normH1_4   = normH1_4  + (deriv1_val-ref)**2*h1
end do

call sll_delete(ad1d)

!print *, '***********************************************************'
!print *, '              Neumann-Dirichlet case '
!print *, '***********************************************************'

do i=1,NPTS1
   eta1         = eta1_pos(i)
   x(i)         = cos(2.0_f64*sll_pi*eta1)
   reference(i) = cos(2.0_f64*sll_pi*eta1)
end do
xprime(1) = 0.0
xprime(2) = 0.0
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_NEUMANN,SLL_DIRICHLET,SPL_DEG)

call set_values_at_boundary1d(ad1d,                 &
                              value_left=1.0_f64,   &
                              value_right=1.0_f64,  &
                              slope_left=xprime(1), &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc4 = 0.0_f64
acc4_der1 = 0.0_f64
normL2_4 = 0.0_f64
normH1_4 = 0.0_f64
do i=1,NPTS1
   eta1       = eta1_pos(i)
   node_val   = ad1d%interpolate_value(eta1)
   ref        = cos(2.0_f64*sll_pi*eta1)
   acc4       = acc4 + abs(node_val-ref)
   normL2_4   = normL2_4  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
   acc4_der1  = acc4_der1 + abs(deriv1_val-ref)
   normH1_4   = normH1_4  + (deriv1_val-ref)**2*h1
end do

call sll_delete(ad1d)

!print *, '***********************************************************'
!print *, '              Neumann-Neumann case '
!print *, '***********************************************************'

do i=1,NPTS1
   eta1         = eta1_pos(i)
   x(i)         = cos(2.0_f64*sll_pi*eta1)
   reference(i) = cos(2.0_f64*sll_pi*eta1)
end do
xprime(1) = 0.0
xprime(2) = 0.0
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_NEUMANN,SLL_NEUMANN,SPL_DEG)

call set_values_at_boundary1d(ad1d,                 &
                              value_left=1.0_f64,   &
                              value_right=1.0_f64,  &
                              slope_left=xprime(1), &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc5 = 0.0_f64
acc5_der1 = 0.0_f64
normL2_5 = 0.0_f64
normH1_5 = 0.0_f64
do i=1,NPTS1
   eta1       = eta1_pos(i)
   node_val   = ad1d%interpolate_value(eta1)
   ref        = cos(2.0_f64*sll_pi*eta1)
   acc5       = acc5 + abs(node_val-ref)
   normL2_5   = normL2_5  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
   acc5_der1  = acc5_der1 + abs(deriv1_val-ref)
   normH1_5   = normH1_5  + (deriv1_val-ref)**2*h1
end do

call sll_delete(ad1d)

!print *, '***********************************************************'
!print *, '              Hermite-Neumann case '
!print *, '***********************************************************'

do i=1,NPTS1
   eta1         = eta1_pos(i)
   x(i)         = cos(2.0_f64*sll_pi*eta1)
   reference(i) = cos(2.0_f64*sll_pi*eta1)
end do
xprime(1) = 0.0
xprime(2) = 0.0
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_HERMITE,SLL_NEUMANN,SPL_DEG)

call set_values_at_boundary1d(ad1d,                  &
                              value_left=1.0_f64,    &
                              value_right=1.0_f64,   &
                              slope_left=xprime(1),  &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc6 = 0.0_f64
acc6_der1 = 0.0_f64
normL2_6 = 0.0_f64
normH1_6 = 0.0_f64
do i=1,NPTS1
   eta1       = eta1_pos(i)
   node_val   = ad1d%interpolate_value(eta1)
   ref        = cos(2.0_f64*sll_pi*eta1)
   acc6       = acc6 + abs(node_val-ref)
   normL2_6   = normL2_6  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
   acc6_der1  = acc6_der1 + abs(deriv1_val-ref)
   normH1_6   = normH1_6  + (deriv1_val-ref)**2*h1
end do

call sll_delete(ad1d)

!print *, '***********************************************************'
!print *, '              Neumann-Hermite case '
!print *, '***********************************************************'

do i=1,NPTS1
   eta1         = eta1_pos(i)
   x(i)         = cos(2.0_f64*sll_pi*eta1)
   reference(i) = cos(2.0_f64*sll_pi*eta1)
end do
xprime(1) = 0.0
xprime(2) = 0.0
eta1_prime(1) = eta1_pos(1)
eta1_prime(2) = eta1_pos(NPTS1)

call ad1d%initialize(NPTS1,X1MIN,X1MAX,SLL_NEUMANN,SLL_HERMITE,SPL_DEG)

call set_values_at_boundary1d(ad1d,                 &
                              value_left=1.0_f64,   &
                              value_right=1.0_f64,  &
                              slope_left=xprime(1), &
                              slope_right=xprime(2))

call ad1d%compute_interpolants(x)

acc7 = 0.0_f64
acc7_der1 = 0.0_f64
normL2_7 = 0.0_f64
normH1_7 = 0.0_f64
do i=1,NPTS1
   eta1       = eta1_pos(i)
   node_val   = ad1d%interpolate_value(eta1)
   ref        = cos(2.0_f64*sll_pi*eta1)
   acc7        = acc7 + abs(node_val-ref)
   normL2_7    = normL2_7  + (node_val-ref)**2*h1
   deriv1_val = ad1d%interpolate_derivative_eta1(eta1)   
   ref        = -2.0_f64*sll_pi*sin(2.0_f64*sll_pi*eta1)
   acc7_der1  = acc7_der1 + abs(deriv1_val-ref)
   normH1_7   = normH1_7  + (deriv1_val-ref)**2*h1
end do

call sll_delete(ad1d)


print*,'--------------------------------------------'
print*,' Average error in nodes'
print*,'--------------------------------------------'
print*,'periodic = ', acc/(NPTS1)
print*,'dirichlet = ', acc1/(NPTS1)
print*,'dirichlet non homogene = ', acc2/(NPTS1)
print*,'Hermite-Hermite = ', acc3/(NPTS1)
print*,'Neumann-Dirichlet = ', acc4/(NPTS1)
print*,'Neumann-Neumann = ', acc5/(NPTS1)
print*,'Hermite-Neumann = ', acc6/(NPTS1)
print*,'Neumann-Hermite = ', acc7/(NPTS1)
print*,' '
print*,'--------------------------------------------'
print*,' Average error in nodes first derivative'
print*,'--------------------------------------------'
print*,'periodic=',acc_der1/(NPTS1)
print*,'dirichlet=',acc1_der1/(NPTS1)
print*,'dirichlet non homogene=',acc2_der1/(NPTS1)
print*,'Hermite-Hermite=',acc3_der1/(NPTS1)
print*,'Neumann-Dirichlet=',acc4_der1/(NPTS1)
print*,'Neumann-Neumann=',acc5_der1/(NPTS1)
print*,'Hermite-Neumann=',acc6_der1/(NPTS1)
print*,'Neumann-Hermite=',acc7_der1/(NPTS1)
print*, ' '
print*,'--------------------------------------------'
print*,' Norm L2 error '
print*,'--------------------------------------------'
print*,'periodic =', sqrt(normL2_0), h1**(SPL_DEG)
print*,'dirichlet =', sqrt(normL2_1),h1**(SPL_DEG)
print*,'dirichlet non homogene =', sqrt(normL2_2),h1**(SPL_DEG)
print*,'Hermite-Hermite =', sqrt(normL2_3),h1**(SPL_DEG)
print*,'Neumann-dirichlet =', sqrt(normL2_4),h1**(SPL_DEG)
print*,'Neumann-Neumann =', sqrt(normL2_5),h1**(SPL_DEG)
print*,'Hermite-Neumann =', sqrt(normL2_6),h1**(SPL_DEG)
print*,'Neumann-Hermite =', sqrt(normL2_7),h1**(SPL_DEG)
print*, ' '
print*,'--------------------------------------------'
print*,' Norm H1 error '
print*,'--------------------------------------------'
print*,'periodic  =', sqrt(normH1_0),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
print*,'dirichlet =', sqrt(normH1_1),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
print*,'dirichlet non homogene =', sqrt(normH1_2),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
print*,'Hermite-Hermite =',sqrt(normH1_3),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
print*,'Neumann-Dirichlet =',sqrt(normH1_4),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
print*,'Neumann-Neumann =',sqrt(normH1_5),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
print*,'Hermite-Neumann =',sqrt(normH1_6),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2
print*,'Neumann-Hermite =',sqrt(normH1_7),h1**(SPL_DEG-2)*(2.0_f64*sll_pi)**2

if(( sqrt(normL2_0) <= h1**(SPL_DEG+1)) .AND. &
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
deallocate(eta1_pos)
deallocate(eta1_prime)

end program unit_test
