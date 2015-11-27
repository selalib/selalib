module woodbury
implicit none

real(8) , allocatable :: u(:,:) 
real(8) , allocatable :: v(:,:)
real(8) , allocatable :: h(:,:)

contains

subroutine woodbury_fac(n, k, q)
implicit none

integer :: n
integer :: k

real(8) :: q(2*k+1,n)
integer :: i
integer :: j
integer :: l
integer :: kp1

real(8) :: g(k)

real(8) :: work(k*k)
integer :: iflag
integer :: info
integer :: jpiv(k)

kp1 = k+1

allocate(u(n,k))
allocate(v(n,k))
allocate(h(k,k))

u = 0.0_8
v = 0.0_8
do l = 1, k
  g(l)   = -q(kp1,l)  ! Arbitrary gamma is set using diagonal term of A
  u(l,l) = g(l)
  v(l,l) = 1.0_8 
end do

do j = 1, k
  l = 0
  do i = n-k+j,n
    l = l+1
    u(i,j) = q(l,j)
    v(i,j) = q(kp1+kp1-l,i)/g(j)
  end do
end do

do j = 1, k
  q(kp1,j) = q(kp1,j) - g(j) 
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


end subroutine woodbury_fac

subroutine woodbury_slv(n, k, q, b, x)
implicit none
integer :: n
integer :: k
real(8) :: x(n), b(n)
real(8) :: q(2*k+1,n)
integer :: kp1

kp1 = k+1
!Solve A.y = b
x = b
call banslv ( q, k+kp1, n, k, k, x )
x = x - matmul(u,matmul(h,matmul(transpose(v),x)))

end subroutine woodbury_slv

end module woodbury
