module woodbury
implicit none

type :: woodbury_solver
  real(8) , allocatable :: u(:,:) 
  real(8) , allocatable :: v(:,:)
  real(8) , allocatable :: h(:,:)
end type woodbury_solver

contains

subroutine woodbury_fac(s, n, k, q)
type(woodbury_solver) :: s
integer :: n
integer :: k

real(8) :: q(2*k+1,n)
integer :: i
integer :: j
integer :: l
integer :: kp1

real(8) :: g(k)

integer :: iflag
real(8), allocatable :: qq(:,:)

kp1 = k+1

allocate(s%u(n,k))
allocate(s%v(n,k))
allocate(s%h(k,k))

s%u = 0.0_8
s%v = 0.0_8
do l = 1, k
  g(l)   = -q(kp1,l)  ! Arbitrary gamma is set using diagonal term of A
  s%u(l,l) = g(l)
  s%v(l,l) = 1.0_8 
end do

do j = 1, k
  l = 0
  do i = n-k+j,n
    l = l+1
    s%u(i,j) = q(l,j)
    s%v(i,j) = q(kp1+kp1-l,i)/g(j)
  end do
end do

do j = 1, k
  q(kp1,j) = q(kp1,j) - g(j) 
end do
do j = n-k+1,n
  l = n-j 
  do i = n-k+1,n
    l = l+1
    q(l,j) = q(l,j) - sum(s%u(i,:)*s%v(j,:)) !This sum can be optimized
  end do
end do

!Factorize the matrix A
call banfac ( q, k+kp1, n, k, k, iflag )

!Solve A.z = u
do l = 1, k
  call banslv ( q, k+kp1, n, k, k, s%u(:,l) )
end do

!compute the matrix H = inverse(1+t(v).z)
s%h = 0.0_8
do i = 1, k
  s%h(i,i) = 1.0_8
end do
s%h = s%h + matmul(transpose(s%v),s%u)

!Inverse H
allocate(qq(k+k-1,k)); qq = 0.0_8
do j = 1, k
  l = 0
  do i = -k+1,k-1
    l=l+1
    if(i+j >=1 .and. i+j<=k) qq(l,j) = s%h(i+j,j)
  end do
end do

call banfac ( qq, k+k-1, k, k-1, k-1, iflag )
do j = 1, k
  do i = 1, k
    if ( i==j) then
      s%h(i,j) = 1.0_8
    else
      s%h(i,j) = 0.0_8
    end if
  end do
  call banslv ( qq, k+k-1, k, k-1, k-1, s%h(:,j) )
end do

end subroutine woodbury_fac

subroutine woodbury_slv(s, n, k, q, x)
type(woodbury_solver) :: s

integer :: n
integer :: k
real(8) :: x(n)
real(8) :: q(2*k+1,n)
integer :: kp1

kp1 = k+1
call banslv ( q, k+kp1, n, k, k, x )
x = x - matmul(s%u,matmul(s%h,matmul(transpose(s%v),x)))

end subroutine woodbury_slv

end module woodbury
