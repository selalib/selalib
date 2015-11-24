program test_cyclic_banded_solver
implicit none

integer, parameter :: n = 10
integer, parameter :: k = 2
real(8) :: x(n)
real(8) :: y(n)
real(8) :: b(n)
real(8) :: q(2*k+1,n)
real(8) :: a(n,n)
real(8) :: u(n,k)
real(8) :: v(n,k)
real(8) :: z(n,k)
integer :: iflag
integer :: i
integer :: j, l
integer :: kp1
real(8) :: fact

real(8) :: m(n,n)
real(8) :: g(k)

kp1 = k+1
do i = 1, n
  b(i) = 1.0_8 
end do

m = 0.0_8
call random_number(a)
write(*,*) "m="
do j = 1, n
  do i = -k,k
    l=modulo(j+i-1,n)+1 
    m(l,j) = a(l,j)
  end do
  write(*,"(11f7.3)") m(:,j)
end do
write(*,*)

a = 0.0
do j = 1, n
  do i = -k,k
    if (i+j>=1 .and. i+j<=n) a(i+j,j) = m(i+j,j)
  end do
end do
print*, "a="
call print_matrix(a)

!set u and v vectors and modify a
u = 0.0_8
v = 0.0_8
do l = 1, k
  g(l)   = -m(l,l)
  u(l,l) = g(l)
  v(l,l) = 1.0 
end do

do j = 1, k
  a(j,j) = a(j,j) - g(j) 
  do i = n-j,n
    u(i,j) = m(i,j)
    v(i,j) = m(j,i)/g(j)
  end do
end do

!a(n-1,n-1) = a(n-1,n-1) - u(n-1,1)*v(n-1,1)
!a(n-1,n  ) = a(n-1,n  ) - u(n-1,1)*v(n  ,1)
!a(n  ,n-1) = a(n  ,n-1) - u(n  ,1)*v(n-1,1)
!a(n  ,n  ) = a(n  ,n  ) - u(n  ,1)*v(n  ,1) - u(n,2)*v(n,2)

do j = n-k+1,n
  do i = n-k+1,n
    a(i,j) = a(i,j) - sum(u(i,:)*v(j,:))
  end do
end do

print*, 'u='
call print_matrix(u)
print*, 'v='
call print_matrix(v)
print*, 'w='
call print_matrix(m - (a + matmul(u,transpose(v))))
stop


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



!do l = 1, k
!
!  q(k+1,1) = a(1,1) - u(1,l) 
!  q(k+1,n) = a(n,n) - u(n,l)*v(n,l)
!
!  call banfac ( q, k+kp1, n, k, k, iflag )
!  y = b
!  call banslv ( q, k+kp1, n, k, k, y )
!
!  z = u
!  call banslv ( q, k+kp1, n, k, k, z )
!
!end do

!fact = (dot_product(y,v))/(1.0_8+dot_product(v,z))
!x = y - fact * z 
!
!write(*,"(' x = ', g15.3)") sum(b - matmul(a,x))


contains


subroutine print_matrix(m)
real(8), intent(in)  :: m(:,:)
integer :: j
character(len=20) :: display_format

write(display_format, "('(''|''',i2,a,''' |'')')")size(m,1),'f7.3'
do j = 1, size(m,2)
  write(*,display_format) m(:,j)
end do
write(*,*)
end subroutine print_matrix

end program test_cyclic_banded_solver
