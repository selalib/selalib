module sll_m_fcisl_toroidal
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_common_array_initializers, only: &
    sll_scalar_initializer_2d

  use sll_m_constants, only: &
    sll_pi

  use sll_m_hermite_interpolation_2d, only: &
    compute_hermite_derivatives_periodic1, &
    compute_w_hermite, &
    localize_per

  use sll_m_lagrange_interpolation, only: &
    lagrange_interpolate

  implicit none

  public :: &
    compute_analytic_field, &
    compute_euler_field, &
    compute_feet_analytic, &
    compute_feet_euler, &
    compute_feet_rk4, &
    compute_inverse_invr_integral, &
    compute_invr_integral, &
    compute_linspace, &
    compute_modulo_vect, &
    compute_modulo_vect2d_inplace, &
    compute_rk4_field, &
    compute_time_points, &
    interpolate2d_toroidal

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  !sll_real64 :: ph
  !sll_real64 :: th
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
  procedure(sll_scalar_initializer_2d), pointer :: func1
  sll_real64 :: params1(2)
  procedure(sll_scalar_initializer_2d), pointer :: func2
  sll_real64 :: params2(3)
  sll_real64 :: dt  
  !sll_real64 :: bigr
  !sll_real64 :: ph
  !sll_real64 :: th
  sll_int32 :: i
  sll_real64 :: k1(4)
  sll_real64 :: k2(4)
  
  func1 => toroidal_1
  params1(1) = smallr
  params1(2) = psipr
  func2 => toroidal_2
  params2(1) = smallr
  params2(2) = R0
  params2(3) = F0 
  theta(1) = theta0
  phi(1) = phi0
  do i=1,num_time_points-1
    dt = time_points(i+1)-time_points(i)
    k1(1) = func1(theta(i),phi(i),params1)
    k2(1) = func2(theta(i),phi(i),params2)
    k1(2) = func1(theta(i)+0.5_f64*dt*k1(1),phi(i)+0.5_f64*dt*k2(1),params1)
    k2(2) = func2(theta(i)+0.5_f64*dt*k1(1),phi(i)+0.5_f64*dt*k2(1),params2)
    k1(3) = func1(theta(i)+0.5_f64*dt*k1(2),phi(i)+0.5_f64*dt*k2(2),params1)
    k2(3) = func2(theta(i)+0.5_f64*dt*k1(2),phi(i)+0.5_f64*dt*k2(2),params2)
    k1(4) = func1(theta(i)+dt*k1(3),phi(i)+dt*k2(3),params1)
    k2(4) = func2(theta(i)+dt*k1(3),phi(i)+dt*k2(3),params2)   
    theta(i+1) = theta(i)+(dt/6._f64)*(k1(1)+2._f64*(k1(2)+k1(3))+k1(4))
    phi(i+1) = phi(i)+(dt/6._f64)*(k2(1)+2._f64*(k2(2)+k2(3))+k2(4))
  enddo

end subroutine compute_rk4_field 

function toroidal_1(theta,phi,params) result(res)
  sll_real64 :: res
  sll_real64, intent(in) :: theta
  sll_real64, intent(in) :: phi
  sll_real64, dimension(:), intent(in), optional :: params

  sll_real64 :: psipr
  sll_real64 :: smallr
  smallr = params(1)
  psipr = params(2)
  res = -psipr/smallr  
  return !PN ADD TO PREVENT WARNING
  print*, phi, theta
end function toroidal_1

function toroidal_2(theta,phi,params) result(res)
  sll_real64 :: res
  sll_real64, intent(in) :: theta
  sll_real64, intent(in) :: phi
  sll_real64, dimension(:), intent(in), optional :: params

  sll_real64 :: bigr
  sll_real64 :: F0
  sll_real64 :: R0
  sll_real64 :: smallr
  
  smallr = params(1)
  R0 = params(2)
  F0 = params(3)
  bigr = R0+smallr*cos(theta)
  res=F0/bigr

  return !PN ADD TO PREVENT WARNING
  print*, phi
end function toroidal_2



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
  !sll_real64 :: dt
  
  !sll_real64 :: bigr
  !sll_real64 :: ph
  !sll_real64 :: th
  sll_int32 :: i

  
  do i=1,num_time_points
    theta(i) = theta0-psipr/smallr*(time_points(i)-time_points(1))
    phi(i) = phi0 -F0*smallr/psipr*(compute_invR_integral(R0,smallr,theta(i))-compute_invR_integral(R0,smallr,theta0))
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
  k = real(floor(xx),f64)
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


function compute_inverse_invR_integral(R0,smallr,x) result(res)
  sll_real64 :: res
  sll_real64, intent(in) :: R0
  sll_real64, intent(in) :: smallr
  sll_real64, intent(in) :: x
  
  sll_real64 :: xx
  sll_real64 :: k
  sll_real64 :: alpha
  !computes inverse of function int(1/(R0+smallr*cos(u)),u=0..x)
  !we suppose that R0>r
  
  xx = 0.5_f64*x*sqrt(R0**2-smallr**2)
  xx = 0.5_f64+xx/sll_pi
  k = real(floor(xx),f64)
  alpha = xx-k
  if(alpha<1.e-12)then 
    res = -0.5_f64*sll_pi
  else
    res = tan((-0.5_f64+alpha)*sll_pi)
    res = atan(sqrt(R0**2-smallr**2)/(R0-smallr)*res)
  endif
  res = res+k*sll_pi  
  res = 2._f64*res
  
end function compute_inverse_invR_integral



subroutine compute_modulo(x,L)
  sll_real64, intent(inout) :: x
  sll_real64, intent(in) :: L
  
  x = x/L+1._f64
  x = x-real(floor(x),f64)
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


subroutine compute_modulo_vect2d_inplace(inout,num_points1,num_points2,L)
  sll_real64, dimension(:,:), intent(inout) :: inout
  sll_int32, intent(in) :: num_points1
  sll_int32, intent(in) :: num_points2
  sll_real64, intent(in) :: L
  
  sll_int32 :: i1
  sll_int32 :: i2
  
  if(size(inout,1)<num_points1)then
    print *,'#problem of size1 of inout',size(inout,1),num_points1
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif
  if(size(inout,2)<num_points2)then
    print *,'#problem of size2 of inout',size(inout,2),num_points2
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif
  
  do i2=1, num_points2
    do i1=1,num_points1
      call compute_modulo(inout(i1,i2),L)
    enddo  
  enddo
  
end subroutine compute_modulo_vect2d_inplace



subroutine compute_linspace(array,xmin,xmax,Npts)
  sll_real64, dimension(:), intent(out) :: array
  sll_real64, intent(in) :: xmin
  sll_real64, intent(in) :: xmax
  sll_int32, intent(in) :: Npts
  
  sll_int32 :: i
  sll_real64 :: dx
  
  dx = (xmax-xmin)/real(Npts-1,f64)
  
  do i=1,Npts
    array(i) = xmin+real(i-1,f64)*dx
  enddo
  
end subroutine compute_linspace

subroutine compute_feet_euler(theta_in,phi_in,num_points1,num_points2,dt,params,theta_out,phi_out)
  sll_real64, dimension(:), intent(in) :: theta_in
  sll_real64, dimension(:), intent(in) :: phi_in
  sll_int32, intent(in) :: num_points1
  sll_int32, intent(in) :: num_points2
  sll_real64, intent(in) :: dt
  sll_real64, dimension(:), intent(in) :: params
  sll_real64, dimension(:,:), intent(out) :: theta_out
  sll_real64, dimension(:,:), intent(out) :: phi_out

  sll_real64 :: R0
  sll_real64 :: psipr
  sll_real64 :: F0
  sll_real64 :: smallr

  sll_real64 :: bigr
  !sll_real64 :: ph
  !sll_real64 :: th
  sll_int32 :: i
  sll_int32 :: j



  if(size(params)<4)then
    print *,'#size of params not good',size(params)
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif

  R0 = params(1)
  psipr = params(2)
  F0 = params(3)
  smallr = params(4)
    
  do j=1,num_points2
    do i=1,num_points1
      bigr = R0+smallr*cos(theta_in(i))
      theta_out(i,j) = theta_in(i)-psipr/smallr*dt
      phi_out(i,j) = phi_in(j)+F0/bigr*dt
    enddo
  enddo  

end subroutine compute_feet_euler


subroutine compute_feet_rk4(theta_in,phi_in,num_points1,num_points2,dt,params,theta_out,phi_out)
  sll_real64, dimension(:), intent(in) :: theta_in
  sll_real64, dimension(:), intent(in) :: phi_in
  sll_int32, intent(in) :: num_points1
  sll_int32, intent(in) :: num_points2
  sll_real64, intent(in) :: dt
  sll_real64, dimension(:), intent(in) :: params
  sll_real64, dimension(:,:), intent(out) :: theta_out
  sll_real64, dimension(:,:), intent(out) :: phi_out

  sll_real64 :: R0
  sll_real64 :: psipr
  sll_real64 :: F0
  sll_real64 :: smallr

  !sll_real64 :: bigr
  !sll_real64 :: ph
  !sll_real64 :: th
  sll_int32 :: i
  sll_int32 :: j

  procedure(sll_scalar_initializer_2d), pointer :: func1
  sll_real64 :: params1(2)
  procedure(sll_scalar_initializer_2d), pointer :: func2
  sll_real64 :: params2(3)
  sll_real64 :: k1(4)
  sll_real64 :: k2(4)




  if(size(params)<4)then
    print *,'#size of params not good',size(params)
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif

  R0 = params(1)
  psipr = params(2)
  F0 = params(3)
  smallr = params(4)

  
  func1 => toroidal_1
  params1(1) = smallr
  params1(2) = psipr
  func2 => toroidal_2
  params2(1) = smallr
  params2(2) = R0
  params2(3) = F0 
    
    
  do j=1,num_points2
    do i=1,num_points1
      k1(1) = func1(theta_in(i),phi_in(j),params1)
      k2(1) = func2(theta_in(i),phi_in(j),params2)
      k1(2) = func1(theta_in(i)+0.5_f64*dt*k1(1),phi_in(j)+0.5_f64*dt*k2(1),params1)
      k2(2) = func2(theta_in(i)+0.5_f64*dt*k1(1),phi_in(j)+0.5_f64*dt*k2(1),params2)
      k1(3) = func1(theta_in(i)+0.5_f64*dt*k1(2),phi_in(j)+0.5_f64*dt*k2(2),params1)
      k2(3) = func2(theta_in(i)+0.5_f64*dt*k1(2),phi_in(j)+0.5_f64*dt*k2(2),params2)
      k1(4) = func1(theta_in(i)+dt*k1(3),phi_in(j)+dt*k2(3),params1)
      k2(4) = func2(theta_in(i)+dt*k1(3),phi_in(j)+dt*k2(3),params2)   
      theta_out(i,j) = theta_in(i)+(dt/6._f64)*(k1(1)+2._f64*(k1(2)+k1(3))+k1(4))
      phi_out(i,j) = phi_in(j)+(dt/6._f64)*(k2(1)+2._f64*(k2(2)+k2(3))+k2(4))
    enddo
  enddo  

end subroutine compute_feet_rk4



subroutine compute_feet_analytic(theta_in,phi_in,num_points1,num_points2,dt,params,theta_out,phi_out)
  sll_real64, dimension(:), intent(in) :: theta_in
  sll_real64, dimension(:), intent(in) :: phi_in
  sll_int32, intent(in) :: num_points1
  sll_int32, intent(in) :: num_points2
  sll_real64, intent(in) :: dt
  sll_real64, dimension(:), intent(in) :: params
  sll_real64, dimension(:,:), intent(out) :: theta_out
  sll_real64, dimension(:,:), intent(out) :: phi_out

  sll_real64 :: R0
  sll_real64 :: psipr
  sll_real64 :: F0
  sll_real64 :: smallr

  sll_real64 :: bigr
  !sll_real64 :: ph
  !sll_real64 :: th
  sll_int32 :: i
  sll_int32 :: j



  if(size(params)<4)then
    print *,'#size of params not good',size(params)
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif

  R0 = params(1)
  psipr = params(2)
  F0 = params(3)
  smallr = params(4)
    
  do j=1,num_points2
    do i=1,num_points1
      bigr = R0+smallr*cos(theta_in(i))
      theta_out(i,j) = theta_in(i)-psipr/smallr*dt
      phi_out(i,j) = phi_in(j) &
        -F0*smallr/psipr*( &
        compute_invR_integral(R0,smallr,theta_out(i,j)) &
        -compute_invR_integral(R0,smallr,theta_in(i)))
    enddo
  enddo  

end subroutine compute_feet_analytic

function interpolate2d_toroidal( &
  num_points1, &
  num_points2, &
  data_in, &
  eta1, &
  eta2, &
  params) &
  result(data_out)
  sll_int32, intent(in) :: num_points1
  sll_int32, intent(in) :: num_points2
  sll_real64, dimension(:,:), intent(in) :: eta1
  sll_real64, dimension(:,:), intent(in) :: eta2
  sll_real64, dimension(:,:), intent(in) :: data_in
  sll_real64, dimension(num_points1,num_points2) :: data_out
  sll_real64, dimension(:), intent(in) :: params
  
  sll_real64 :: R0
  sll_real64 :: psipr
  sll_real64 :: F0
  sll_real64 :: smallr
  sll_int32 :: hermite_p
  sll_int32 :: lag_r
  sll_int32 :: lag_s
  sll_int32 :: lag_p
  sll_int32 :: ierr
  sll_real64, dimension(:,:,:), allocatable :: hermite_buf
  sll_int32 :: i
  sll_int32 :: j
  sll_real64 :: eta1_min
  sll_real64 :: eta1_max
  sll_real64 :: eta2_min
  sll_real64 :: eta2_max
  !sll_real64 :: eta1_star
  !sll_real64 :: eta2_star
  sll_real64, dimension(:), allocatable :: lag_buf
  sll_real64, dimension(:), allocatable :: lag_x
  sll_real64 :: eta1_loc 
  sll_real64 :: eta2_loc 
  sll_int32 :: i0
  sll_int32 :: j0
  sll_int32 :: j0_loc
  sll_int32 :: jj
  sll_real64 :: w(4)
  sll_real64 :: tmp1
  sll_real64 :: tmp2
  !we define aligned interpolation
  !Lagrange interpolation in aligned direction (lag_r,lag_s)
  !Hermite interpolation in theta direction (hermite_p)
  !alignement with respect to magnetic field lines
  !given analytically 
  
  R0 = params(1)
  psipr = params(2)
  F0 = params(3)
  smallr = params(4)
  hermite_p = int(params(5),i32)
  lag_r = int(params(6),i32)
  lag_s = int(params(7),i32)
  eta1_min = params(8)
  eta1_max = params(9)
  eta2_min = params(10)
  eta2_max = params(11)
  
  lag_p = lag_s-lag_r
  
  SLL_ALLOCATE(hermite_buf(3,num_points1,num_points2),ierr)
  call compute_hermite_derivatives_periodic1( &
    data_in, &
    num_points1, &
    num_points2, &
    hermite_p, &
    hermite_buf)
  
  !print *,'hermite_buf=',maxval(abs(hermite_buf(1,:,:)-1._f64)), &
  !  maxval(abs(hermite_buf(2,:,:))),maxval(abs(hermite_buf(3,:,:)))
  !stop
  
  SLL_ALLOCATE(lag_buf(lag_r:lag_s),ierr)
  SLL_ALLOCATE(lag_x(lag_r:lag_s),ierr)
  do i=lag_r,lag_s
    lag_x(i) = real(i,f64)
  enddo

  tmp1 = (eta2_max-eta2_min)/real(num_points2-1,f64)*(-psipr)/(F0*smallr)
  do j=1,num_points2
    do i=1,num_points1
      eta2_loc = eta2(i,j)
      call localize_per(j0,eta2_loc,eta2_min,eta2_max,num_points2-1)
      tmp2 = compute_invR_integral(R0,smallr,eta1(i,j))
      tmp2 = tmp2-eta2_loc*tmp1
      do jj=lag_r,lag_s
        j0_loc = modulo(j0+jj+num_points2-1,num_points2-1)+1
        eta1_loc = compute_inverse_invR_integral(R0,smallr,real(jj,f64)*tmp1+tmp2)
        !phi(t) = phi(0)-F0r/psipr*(L(th(t))-L(th(0)))
        !Hermite interpolation
        call localize_per(i0,eta1_loc,eta1_min,eta1_max,num_points1-1)
        i0=i0+1
        w(1)=(2._f64*eta1_loc+1)*(1._f64-eta1_loc)*(1._f64-eta1_loc)
        w(2)=eta1_loc*eta1_loc*(3._f64-2._f64*eta1_loc)
        w(3)=eta1_loc*(1._f64-eta1_loc)*(1._f64-eta1_loc)
        w(4)=eta1_loc*eta1_loc*(eta1_loc-1._f64)
        
        !if(abs(eta2_loc)<1.e-12)then
        !  print *,eta1_loc,eta2_loc,j0,jj
        !endif
        lag_buf(jj) = w(1)*hermite_buf(1,i0,j0_loc) &
          +w(2)*hermite_buf(1,i0+1,j0_loc) &
          +w(3)*hermite_buf(2,i0,j0_loc) &
          +w(4)*hermite_buf(3,i0,j0_loc)
        !if(abs(lag_buf(jj)-1._f64)>1.e-12)then  
        !  print *,i0,j0_loc,lag_buf(jj)-1._f64,hermite_buf(1,i0,j0_loc),hermite_buf(2,i0,j0_loc),hermite_buf(3,i0,j0_loc),&
        !  hermite_buf(1,i0+1,j0_loc)
        !  stop
        !endif  
      enddo
      
      data_out(i,j) = lagrange_interpolate( &
        eta2_loc, &
        lag_p, &
        lag_x(lag_r:lag_s), &
        lag_buf(lag_r:lag_s))
    enddo
  enddo
end function interpolate2d_toroidal

subroutine compute_w_hermite_aligned( &
  w, &
  w_cell, &
  num_points1, &
  r, &
  s, &
  eta1_pos, &
  eta1_min, &
  eta1_max )
  sll_real64, dimension(:,:,:), intent(out) :: w
  sll_int32, dimension(:,:), intent(out) :: w_cell
  sll_int32, intent(in) :: num_points1
  sll_int32, intent(in) :: r
  sll_int32, intent(in) :: s
  sll_real64, dimension(:,:), intent(in) :: eta1_pos
  sll_real64, intent(in) :: eta1_min
  sll_real64, intent(in) :: eta1_max
  sll_int32 :: i
  !sll_int32 :: j
  sll_int32 :: ell
  sll_real64 :: w_loc_tmp(r:s)
  sll_real64 :: w_loc(s-r)
  sll_int32 :: i0
  sll_real64 :: eta1_loc
  sll_real64 :: w_basis(4)
  sll_int32 :: ii
  character(len=*), parameter :: fun = "compute_w_hermite_aligned" 
  
  
  
  if(r>0 .or. s<0) then
    print *,'#r,s=',r,s
    print *,'#not treated'
    SLL_ERROR(fun,'#bad value of r and s')
  endif
  if(size(w,1)<4)then
    SLL_ERROR(fun,'#bad size1 for w')
  endif
  if(size(w,2)<s-r)then
    SLL_ERROR(fun,'#bad size2 for w')
  endif
  if(size(w,3)<num_points1)then
    SLL_ERROR(fun,'#bad size3 for w')
  endif
  if(size(eta1_pos,1)<s-r)then
    SLL_ERROR(fun,'#bad size1 for eta1_pos')
  endif
  if(size(eta1_pos,2)<num_points1)then
    SLL_ERROR(fun,'#bad size2 for eta1_pos')
  endif
  
  call compute_w_hermite(w_loc_tmp,r,s)
  w_loc(1:-r) = w_loc_tmp(r:-1)
  w_loc(-r+1:s) = w_loc_tmp(1:s)
  
  do i=1,num_points1
    do ell=1,r-s
      eta1_loc = eta1_pos(ell,i)
      call localize_per(i0,eta1_loc,eta1_min,eta1_max,num_points1-1)
      i0=i0+1
      w_cell(ell,i) = i0
      w_basis(1)=(2._f64*eta1_loc+1)*(1._f64-eta1_loc)*(1._f64-eta1_loc)
      w_basis(2)=eta1_loc*eta1_loc*(3._f64-2._f64*eta1_loc)
      w_basis(3)=eta1_loc*(1._f64-eta1_loc)*(1._f64-eta1_loc)
      w_basis(4)=eta1_loc*eta1_loc*(eta1_loc-1._f64)        
      do ii=1,4
        w(ii,ell,i) = w_loc(ell)*w_basis(ii)
      enddo  
    enddo   
  enddo

end subroutine compute_w_hermite_aligned


subroutine compute_hermite_derivatives_aligned( &
  f, &
  num_points1, &
  num_points2, &
  p1, &
  w_aligned_left, &
  w_cell_aligned_left, &
  w_aligned_right, &
  w_cell_aligned_right, &
  buf)
  sll_real64, dimension(:,:), intent(in) :: f
  sll_int32, intent(in) :: num_points1
  sll_int32, intent(in) :: num_points2
  sll_int32, intent(in) :: p1
  sll_real64, dimension(:,:,:), intent(in) :: w_aligned_left
  sll_real64, dimension(:,:,:), intent(in) :: w_aligned_right
  sll_int32, dimension(:,:), intent(in) :: w_cell_aligned_left
  sll_int32, dimension(:,:), intent(in) :: w_cell_aligned_right
  sll_real64, dimension(:,:,:), intent(out) :: buf
  
  sll_real64 :: w_left(-p1/2:(p1+1)/2)
  sll_real64 :: w_right((-p1+1)/2:p1/2+1)
  sll_int32 :: r_left
  sll_int32 :: s_left
  sll_int32 :: r_right
  sll_int32 :: s_right
  sll_int32 :: i
  sll_int32 :: j
  sll_real64 :: tmp
  sll_int32 :: ii
  sll_int32 :: ind
      
  r_left=-p1/2
  s_left=(p1+1)/2
  r_right=(-p1+1)/2
  s_right=p1/2+1   
  call compute_w_hermite(w_left,r_left,s_left)
  if(((2*p1/2)-p1)==0)then
    w_right(r_right:s_right) = w_left(r_left:s_left)
  else
    w_right(r_right:s_right) = -w_left(s_left:r_left:-1)
  endif    
  
  
  !first compute \partial_1f(eta1_i,eta2_j)
  do j=1,num_points2
    do i=1,num_points1
      buf(1,i,j) = f(i,j)  !f(0)
      tmp=0._f64
      do ii=r_left,s_left
        ind=modulo(i+ii-2+num_points1,num_points1-1)+1
        tmp=tmp+w_left(ii)*f(ind,j)
      enddo
      buf(2,i,j)=tmp !df(0)
      tmp=0._f64
      do ii=r_right,s_right
        ind=modulo(i+ii-2+num_points1,num_points1-1)+1
        tmp=tmp+w_right(ii)*f(ind,j)
      enddo
      buf(3,i,j)=tmp !df(1)      
    enddo
  enddo

  !then compute \partial_aligned f(eta1_i,eta2_j)

  return !PN ADD TO PREVENT WARNING
  print*, w_aligned_left,       &
          w_cell_aligned_left,  &
          w_aligned_right,      &
          w_cell_aligned_right
  
end subroutine compute_hermite_derivatives_aligned




end module sll_m_fcisl_toroidal
