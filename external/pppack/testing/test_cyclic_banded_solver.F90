program test_cyclic_banded_solver
implicit none

integer, parameter :: n = 10
integer, parameter :: k = 1
real(8) :: f(n)
real(8) :: gtau(n)
real(8) :: tau(n)
real(8) :: z(n)
real(8) :: q(2*k+1,n)
real(8) :: a(n,n)
real(8) :: u(n)
real(8) :: v(n)

real(8) :: tau_min 
real(8) :: tau_max

integer :: iflag
integer :: i
integer :: j, l
integer :: kp1
real(8) :: fact
real(8) :: alpha, beta, gamma
real(8) :: pi

real(8), allocatable :: df_coeffs(:) 


allocate(df_coeffs(2*k+1))
select case(k)
case(1)
  df_coeffs = [real(8) :: 1.,-2.,1.]
case(2)
  df_coeffs = [real(8) :: -1./12., 4./3.,-5./2.,4./3.,-1./12.]
case default
  stop 'wrong value of k'
end select

pi = 4.0_8 * atan(1.0_8)
tau_min = 0.0_8
tau_max = 2.0_8*pi

kp1 = k+1
do i = 1, n
  tau(i)  = tau_min + (i-1)*(tau_max-tau_min)/(n)
  gtau(i) = sin(tau(i)) * ((tau_max-tau_min)/(n))**2
end do

a = 0.0_8
write(*,*) "a="
do j = 1, n
  do i = -k,k
    a(modulo(j+i-1,n)+1,j) = df_coeffs(i+k+1)
  end do
  write(*,"(11f7.3)") a(:,j)
end do
write(*,*)

q = 0.0_8
do j = 1, n
  l = 0
  do i = -k,k
    l = l+1
    if (i+j >= 1 .and. i+j <= n) q(l,j) = a(i+j,j)
  end do
end do

write(*,*) "q="
do i = 1, 2*k+1
  write(*,"(11f7.3)") q(i,:)
end do

gamma      = - a(1,1) 
alpha      = a(1,n)
beta       = a(n,1)

u    = 0.0_8    
u(1) = gamma
u(n) = alpha

v    = 0.0_8    
v(1) = 1.0_8
v(n) = beta / gamma

q(k+1,1)   = a(1,1) - u(1) 
q(k+1,n)   = a(n,n) - u(n)*v(n)

call banfac ( q, k+kp1, n, k, k, iflag )
f = gtau
call banslv ( q, k+kp1, n, k, k, f )

z = u
call banslv ( q, k+kp1, n, k, k, z )

write(*,"(10f7.3)") z
write(*,"(10f7.3)") v

fact = (f(1)+f(n)*v(n))/(1.0_8+v(1)*z(1)+v(n)*z(n))
!fact = dot_product(v,f)/(1.0_8+dot_product(v,z))
f = f - fact * z 

write(*,"(' f = ', 11f7.3)") f

do i = 1, n
   write(17,*) tau(i), f(i) - sum(f)/real(n,8)
end do

end program test_cyclic_banded_solver
