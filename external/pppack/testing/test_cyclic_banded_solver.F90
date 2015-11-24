program test_cyclic_banded_solver
implicit none

integer, parameter :: n = 10
integer, parameter :: k = 3

real(8) :: x(n)
real(8) :: y(n)
real(8) :: b(n)
real(8) :: q(2*k+1,n)
real(8) :: a(n,n)
real(8) :: u(n,k)
real(8) :: v(n,k)
real(8) :: z(n)
integer :: ifla
integer :: i
integer :: j, l
integer :: kp1
real(8) :: fact

real(8) :: m(n,n)
real(8) :: lu(n,n)
real(8) :: g(k)


real(8)              :: h(k,k)
real(8)              :: work(k*k)
integer, allocatable :: ipiv(:)
integer              :: nrhs
integer              :: iflag
integer              :: info


kp1 = k+1
do i = 1, n
  b(i) = 1.0_8 * i 
end do

!Build the complete system with random coefficients
m = 0.0_8
call random_number(a)
!Create a cyclic banded system m
do j = 1, n
  do i = -k,k
    l=modulo(j+i-1,n)+1 
    m(l,j) = a(l,j)
  end do
end do
write(*,*) "M="
call print_matrix(m)

!Create matrix A without corners terms
a = 0.0_8
do j = 1, n
  do i = -k,k
    if (i+j>=1 .and. i+j<=n) a(i+j,j) = m(i+j,j)
  end do
end do

!General Matrix Factorization with Lapack
allocate(ipiv(n)); ipiv=0
lu = m
call dgetrf(n,n,lu,n,ipiv,info)
x = b
nrhs = 1
call dgetrs('n',n,nrhs,lu,n,ipiv,x,n,info)
write(*,"(' Lapack error = ', g15.3)") sum(abs(b-matmul(m,x)))

!set u and v vectors and modify a
u = 0.0_8
v = 0.0_8
do l = 1, k
  g(l)   = -a(l,l)  ! Arbitrary gamma is set by diagonal of A
  u(l,l) = g(l)
  v(l,l) = 1.0_8 
end do

do j = 1, k
  a(j,j) = a(j,j) - g(j) 
  do i = n-j,n
    u(i,j) = m(i,j)
    v(i,j) = m(j,i)/g(j)
  end do
end do

do j = n-k+1,n
  do i = n-k+1,n
    a(i,j) = a(i,j) - sum(u(i,:)*v(j,:)) !This sum can be optimized
  end do
end do

print*, "A:"
call print_matrix(a)
print*, 'U:'
call print_matrix(u)
print*, 'V:'
call print_matrix(v)
print*, 'M - (A + U.t(V)) : Should be zero '
call print_matrix(m - (a + matmul(u,transpose(v))))

!store banded matrix A in q for banfac
q = 0.0_8
do j = 1, n
  l = 0
  do i = -k,k
    l = l+1
    if (i+j >= 1 .and. i+j <= n) q(l,j) = a(i+j,j)
  end do
end do

write(*,*) "Banded matrix A stored in Q:"
call print_matrix(q)

!Factorize the matrix A
call banfac ( q, k+kp1, n, k, k, iflag )

!Solve A.z = u
do l = 1, k
  z = u(:,l)
  call banslv ( q, k+kp1, n, k, k, z )
  u(:,l) = z
end do
print*,'Solve A.Z=U, Z:'
call print_matrix(u)

print*,'(1.+t(V).Z)'
h = 1.0_8+matmul(transpose(v),u)
call print_matrix(h)

!compute the matrix H = inverse(1+t(v).z)

deallocate(ipiv)
allocate(ipiv(k))
call dgetrf(k,k,h,n,ipiv,info)
call dgetri(k,h,k,ipiv,work,k*k,info)
print*,'H . (1.+t(V).Z) : Should be identity'
call print_matrix(matmul(h,1.0_8+matmul(transpose(v),u)))

!NRHS = 1
!call DGETRS('N',NPTS,NRHS,LU,NPTS,IPIV,X,NPTS,INFO)


!fact = (dot_product(y,v))/(1.0_8+dot_product(v,z))
!x = y - fact * z 
!
!write(*,"(' x = ', g15.3)") sum(b - matmul(a,x))


contains


subroutine print_matrix(m)
real(8), intent(in)  :: m(:,:)
integer :: j
character(len=20) :: display_format

write(display_format, "('(''|''',i2,a,''' |'')')")size(m,1),'f8.4'
do j = 1, size(m,2)
  write(*,display_format) m(:,j)
end do
write(*,*)
end subroutine print_matrix

end program test_cyclic_banded_solver
