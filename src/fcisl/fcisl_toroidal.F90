module sll_fcisl_toroidal_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_constants

implicit none

!we consider i this test a case where magnetic filed
!is not straght lines in theta phi
!there is a change
!however everything can be analitically computed
!we write here some routines
!that are useful for making the test
!which is test_fcisl_toroidal

contains 

subroutine compute_time_points(dt,num_time_points,time_points)
  sll_real64, intent(in) :: dt
  sll_int32, intent(in) :: num_time_points
  sll_real64, dimension(:), intent(out) :: time_points
  
  sll_int32 :: i
  
  if(size(time_points)<num_time_points)then
    print *,'#bad size for time_points'
    print *,'#at line/file:',__LINE__,__FILE__
  endif
  
  do i=1,num_time_points
    time_points(i) = real(i-1,f64)*dt
  enddo
  
end subroutine compute_time_points


subroutine compute_euler_field( &
  R0, &
  time_points, &
  num_time_points, &
  psipr, &
  F0, &
  smallr, &
  theta, &
  phi, &
  theta0, &
  phi0)
  sll_real64, intent(in) :: R0
  sll_int32, intent(in) :: num_time_points
  sll_real64, intent(in) :: psipr
  sll_real64, intent(in) :: F0
  sll_real64, intent(in) :: smallr
  sll_real64, dimension(:), intent(in) :: time_points
  sll_real64, dimension(:), intent(out) :: theta
  sll_real64, dimension(:), intent(out) :: phi
  sll_real64, intent(in) :: theta0
  sll_real64, intent(in) :: phi0
  sll_real64 :: dt
  
  sll_real64 :: bigr
  sll_real64 :: ph
  sll_real64 :: th
  sll_int32 :: i
  
  theta(1) = theta0
  phi(1) = phi0
  do i=1,num_time_points-1
    dt = time_points(i+1)-time_points(i)
    bigr = R0+smallr*cos(theta(i))
    theta(i+1) = theta(i)-psipr/smallr*dt
    phi(i+1) = phi(i)+F0/bigr*dt
  enddo

end subroutine compute_euler_field 


subroutine compute_rk4_field( &
  R0, &
  time_points, &
  num_time_points, &
  psipr, &
  F0, &
  smallr, &
  theta, &
  phi, &
  theta0, &
  phi0)
  sll_real64, intent(in) :: R0
  sll_int32, intent(in) :: num_time_points
  sll_real64, intent(in) :: psipr
  sll_real64, intent(in) :: F0
  sll_real64, intent(in) :: smallr
  sll_real64, dimension(:), intent(in) :: time_points
  sll_real64, dimension(:), intent(out) :: theta
  sll_real64, dimension(:), intent(out) :: phi
  sll_real64, intent(in) :: theta0
  sll_real64, intent(in) :: phi0
  sll_real64 :: dt
  
  sll_real64 :: bigr
  sll_real64 :: ph
  sll_real64 :: th
  sll_int32 :: i
  sll_real64 :: k_theta(4)
  sll_real64 :: k_phi(4)
  
  theta(1) = theta0
  phi(1) = phi0
  do i=1,num_time_points-1
    dt = time_points(i+1)-time_points(i)
    bigr = R0+smallr*cos(theta(i))
    theta(i+1) = theta(i)-psipr/smallr*dt
    phi(i+1) = phi(i)+F0/bigr*dt
!    dt = t[i+1]-t[i]
!    k_theta(1) = dp_H(q[i],p[i])
!    k_phi(1) = 
!    kp[0] = -dq_H(q[i],p[i])
!    kq[1] = dp_H(q[i]+0.5*dt*kq[0],p[i]+0.5*dt*kp[0])
!    kp[1] = -dq_H(q[i]+0.5*dt*kq[0],p[i]+0.5*dt*kp[0])
!    kq[2] = dp_H(q[i]+0.5*dt*kq[1],p[i]+0.5*dt*kp[1])
!    kp[2] = -dq_H(q[i]+0.5*dt*kq[1],p[i]+0.5*dt*kp[1])
!    kq[3] = dp_H(q[i]+dt*kq[2],p[i]+dt*kp[2])
!    kp[3] = -dq_H(q[i]+dt*kq[2],p[i]+dt*kp[2])
!    q[i+1] = q[i]+(dt/6.)*(kq[0]+2.*kq[1]+2.*kq[2]+kq[3]) 
!    p[i+1] = p[i]+(dt/6.)*(kp[0]+2.*kp[1]+2.*kp[2]+kp[3]) 
    
    
    
  enddo

end subroutine compute_rk4_field 







subroutine compute_analytic_field( &
  R0, &
  time_points, &
  num_time_points, &
  psipr, &
  F0, &
  smallr, &
  theta, &
  phi, &
  theta0, &
  phi0)
  sll_real64, intent(in) :: R0
  sll_int32, intent(in) :: num_time_points
  sll_real64, intent(in) :: psipr
  sll_real64, intent(in) :: F0
  sll_real64, intent(in) :: smallr
  sll_real64, dimension(:), intent(in) :: time_points
  sll_real64, dimension(:), intent(out) :: theta
  sll_real64, dimension(:), intent(out) :: phi
  sll_real64, intent(in) :: theta0
  sll_real64, intent(in) :: phi0
  sll_real64 :: dt
  
  sll_real64 :: bigr
  sll_real64 :: ph
  sll_real64 :: th
  sll_int32 :: i

  
  do i=1,num_time_points
    theta(i) = theta0-psipr/smallr*(time_points(i)-time_points(1))
    phi(i) = phi0 +F0*smallr/psipr*(compute_invR_integral(R0,smallr,-theta(i))-compute_invR_integral(R0,smallr,-theta0))
  enddo

end subroutine compute_analytic_field 


function compute_invR_integral(R0,smallr,x) result(res)
  sll_real64 :: res
  sll_real64, intent(in) :: R0
  sll_real64, intent(in) :: smallr
  sll_real64, intent(in) :: x
  
  sll_real64 :: xx
  sll_real64 :: k
  sll_real64 :: alpha
  !computes int(1/(R0+smallr*cos(u)),u=0..x)
  !we suppose that R0>r
  
  xx = 0.5_f64+x/(2._f64*sll_pi)
  k = floor(xx)
  alpha = xx-k
  if(alpha<1.e-12)then 
    res = -0.5_f64*sll_pi
  else
    res = tan((-0.5_f64+alpha)*sll_pi)
    res = atan((R0-smallr)/sqrt(R0**2-smallr**2)*res)
  endif
  res = res+k*sll_pi  
  res = 2._f64/(sqrt(R0**2-smallr**2))*res
  
end function compute_invR_integral


subroutine compute_modulo(x,L)
  sll_real64, intent(inout) :: x
  sll_real64, intent(in) :: L
  
  x = x/L+1._f64
  x = x-floor(x)
  x = L*x
  
end subroutine compute_modulo

subroutine compute_modulo_vect(in,out,num_points,L)
  sll_real64, dimension(:), intent(in) :: in
  sll_real64, dimension(:), intent(out) :: out
  sll_int32, intent(in) :: num_points
  sll_real64, intent(in) :: L
  
  sll_int32 :: i
  
  if(size(in)<num_points)then
    print *,'#problem of size of in',size(in),num_points
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif
  if(size(out)<num_points)then
    print *,'#problem of size of out',size(out),num_points
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif
  
  do i=1,num_points
    out(i) = in(i)
    call compute_modulo(out(i),L)
  enddo
  
end subroutine compute_modulo_vect


end module sll_fcisl_toroidal_module