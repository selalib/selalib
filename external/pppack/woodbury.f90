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
subroutine woodbury_solver(m, q, n, k, b)
implicit none
integer, intent(in) :: n
integer, intent(in) :: k

real(8) :: x(n), b(n)
real(8) :: q(2*k+1,n)
real(8) :: a(n,n)
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

!set u and v vectors and modify a
u = 0.0_8
v = 0.0_8
do l = 1, k
  g(l)   = -m(l,l)  ! Arbitrary gamma is set using diagonal term of A
  u(l,l) = g(l)
  v(l,l) = 1.0_8 
end do

do j = 1, k
  do i = n-k+j-1,n
    u(i,j) = m(i,j)
    v(i,j) = m(j,i)/g(j)
  end do
end do

do j = 1, k
  q(kp1,j) = q(kp1,j) - g(j) 
end do

do j = n-k+1,n
  l = k
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

end subroutine woodbury_solver


