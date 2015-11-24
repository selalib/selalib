program test_cyclic_banded_solver
implicit none

integer, parameter :: n = 10
integer, parameter :: k = 1
real(8) :: x(n)
real(8) :: y(n)
real(8) :: b(n)
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

real(8), allocatable  :: df_coeffs(:) 
real(8)               :: lu(n,n)
integer,  allocatable :: ipiv(:)
integer               :: info
integer               :: nrhs

pi = 4.0_8 * atan(1.0_8)
tau_min = 0.0_8
tau_max = 2.0_8*pi

kp1 = k+1
do i = 1, n
  tau(i)  = tau_min + (i-1)*(tau_max-tau_min)/(n)
  b(i) = 1.0_8 !sin(tau(i)) * ((tau_max-tau_min)/(n))**2
end do

lu = 0.0_8
a  = 0.0_8
call random_number(a)
write(*,*) "a="
do j = 1, n
  do i = -k,k
    l=modulo(j+i-1,n)+1 
    lu(l,j) = a(l,j)
  end do
  write(*,"(11f7.3)") lu(:,j)
end do
write(*,*)

a = lu

!General Matrix Factorization
allocate(ipiv(n)); ipiv=0
call dgetrf(n,n,lu,n,ipiv,info)
x = b
nrhs = 1
call dgetrs('n',n,nrhs,lu,n,ipiv,x,n,info)
write(*,"(' x = ', 11f7.3)") b - matmul(a,x)

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

gamma      = -a(1,1) 
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
y = b
call banslv ( q, k+kp1, n, k, k, y )

z = u
call banslv ( q, k+kp1, n, k, k, z )

write(*,"('x:',10f7.3)") x
write(*,"('z:',10f7.3)") z
write(*,"('v:',10f7.3)") v

write(*,"(' z.v = ', f7.3)") z(1)*v(1)+z(n)*v(n)
write(*,"(' z.v = ', f7.3)") sum(z*v)
write(*,"(' z.v = ', f7.3)") dot_product(z,v)

fact = (y(1)+y(n)*v(n))/(1.0_8+v(1)*z(1)+v(n)*z(n))
x = y - fact * z 

write(*,"(' x = ', 11f7.3)") x - matmul(a,b)

do i = 1, n
   write(17,*) tau(i), x(i) !- sum(x)/real(n,8)
end do

end program test_cyclic_banded_solver
