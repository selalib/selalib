! Solve M.x = b
! Compute the A matrix andf U and V vectors where
! A + U.t(V) = M
! To solve (A+U.t(V))X = b. First solve:
! A.z_p = u_p  (for all columns of U)
! Inverse matrix H
! H = inverse(1 + t(V).Z)
! Solve 
! Ay = b
! and the solution is
! x = y - Z.[H.(t(V).y))]
program test_woodbury_solver
implicit none

integer, parameter :: n = 9
integer, parameter :: k = 3

real(8) :: x(n), b(n)
real(8) :: q(2*k+1,n)
real(8) :: u(n,k) 
real(8) :: v(n,k)
integer :: i
integer :: j
integer :: l
integer :: kp1

real(8) :: m(n,n)
real(8) :: g(k)
real(8) :: h(k,k)

real(8)              :: work(k*k)
integer              :: iflag
integer              :: info
integer              :: jpiv(k)


kp1 = k+1
do i = 1, n
  b(i) = 1.0_8 *i
end do

!Build the complete system with random coefficients
!Create a cyclic banded system m
m = 0.0_8
do i = 1, n
  do j = -k,k
    l=modulo(j+i-1,n)+1 
    m(i,l) = real(i*10+l,8)
  end do
end do

call print_matrix(m)

q = 0.0_8
do j = 1, n
  l = 0
  do i = -k,k
    l = l+1
    q(l,j) = m(modulo(i+j-1,n)+1,j)
  end do
end do
call print_matrix(q)


!set u and v vectors and modify a
u = 0.0_8
v = 0.0_8
do l = 1, k
  g(l)   = -q(kp1,l)  ! Arbitrary gamma is set using diagonal term of A
  u(l,l) = g(l)
  v(l,l) = 1.0_8 
end do

do j = 1, k
  do i = n-k+j-1,n
    print*,i*10+j
    u(i,j) = m(i,j)
    v(i,j) = m(j,i)/g(j)
  end do
end do


do j = 1, k
  q(kp1,j) = m(j,j) - g(j) 
end do
do j = n-k+1,n
  l = n-j 
  do i = n-k+1,n
    l = l+1
    q(l,j) = q(l,j) - sum(u(i,:)*v(j,:)) !This sum can be optimized
  end do
end do


!Factorize the matrix A
call banfac ( q, k+kp1, n, k, k, iflag )

!Solve A.y = b
x = b
call banslv ( q, k+kp1, n, k, k, x )
!Solve A.z = u
do l = 1, k
  call banslv ( q, k+kp1, n, k, k, u(:,l) )
end do

!compute the matrix H = inverse(1+t(v).z)
h = 0.0_8
do i = 1, k
  h(i,i) = 1.0_8
end do
h = h + matmul(transpose(v),u)

call dgetrf(k,k,h,k,jpiv,info)
call dgetri(k,h,k,jpiv,work,k*k,info)
x = x - matmul(u,matmul(h,matmul(transpose(v),x)))

write(*,"(' error = ', g15.3)") sum(b - matmul(m,x))

contains

subroutine print_matrix(m)
real(8), intent(in)  :: m(:,:)
character(len=20) :: display_format

write(display_format, "('(''|''',i2,a,''' |'')')")size(m,2),'f9.3'
do i = 1, size(m,1)
  write(*,display_format) m(i,:)
end do
write(*,*)
end subroutine print_matrix

end program test_woodbury_solver


