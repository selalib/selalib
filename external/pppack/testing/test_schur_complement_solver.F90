! To solve M.x = b


!Notre système peut s’écrire sous forme bloc
!(A B)
!(C D)
!où C et B sont des petites matrices.
!En ecrivant le système on à donc à résoudre
!A x1 + B x2 = b1
!C x1 + D x2 = b2
!On peut alors éliminer x1 de la deuxième equation en soustrayant CA^{-1} fois la premiere:
!(D - C A^{-1}B) x2 = b2 -  C A^{-1} b1
!
!On note alors H = D - C A^{-1}B,  c2=  b2 -  C A^{-1} b1
!
!On obtient donc la solution après avoir assemblé H en résolvant
!H x_2 = c2
!A x1 = b1- B x2
!
!d’ou l’algorithme:
!Factorisation:
!Factoriser A (dgbtrf)
!Resoudre: A Y = B  (Y=A^{-1}B), Y est une petite matrice. On peut appeler dgbtrs  avec plusieurs second membres
!Calculer H = D - CY
!Factoriser H
!
!Solve:
!Calculer  c2=b2 -  C A^{-1} b1 (en résolvant  A z2=b1, puis b2=C z2) 
!Resoudre H x_2 = c2 puis A x1 = b1- B x2


program test_cyclic_banded_solver
implicit none

integer, parameter :: n = 9
integer, parameter :: k = 3

real(8) :: x(n)
real(8) :: y(n)
real(8) :: b(n)
real(8) :: q(2*k+1,n)
real(8) :: one(k,k)
real(8) :: a(n,n)
real(8) :: u(n,k) 
real(8) :: v(n,k)
real(8) :: w(n,n)
real(8) :: p(n,n)
real(8) :: z(n,k)
integer :: ifla
integer :: i
integer :: j, l
integer :: kp1
real(8) :: fact

real(8) :: m(n,n)
real(8) :: im(n,n)
real(8) :: lu(n,n)
real(8) :: f(k)
real(8) :: g(k)
real(8) :: h(k,k)

real(8)              :: work(k*k)
integer, allocatable :: ipiv(:)
integer              :: nrhs
integer              :: iflag
integer              :: info
integer              :: ldab
real(8), allocatable :: ab(:,:)
integer              :: jpiv(k)

one = 0.0_8
do i = 1, k
  one(i,i) = 1.0_8
end do

kp1 = k+1
do i = 1, n
  b(i) = 1.0_8 *i
end do

!Build the complete system with random coefficients
call random_number(p) 
!Create a cyclic banded system m
m = 0.0_8
do i = 1, n
  do j = -k,k
    l=modulo(j+i-1,n)+1 
    m(i,l) = p(i,l)
  end do
end do
write(*,*) "M="
call print_matrix(m)


!General Matrix Factorization with Lapack
allocate(ipiv(n)); ipiv=0
lu = m
call dgetrf(n,n,lu,n,ipiv,info)
x = b
nrhs = 1
call dgetrs('n',n,nrhs,lu,n,ipiv,x,n,info)
write(*,"(' Lapack error = ', g15.3)") sum(abs(b-matmul(m,x)))
call print_vector(x)


!Create matrix A without corners terms
a = 0.0_8
do j = 1,n
  do i = max(1,j-k), min(n,j+k)
    a(i,j) = m(i,j)
  end do
end do

!set u and v vectors and modify a
u = 0.0_8
v = 0.0_8
do l = 1, k
  g(l)   = -a(l,l)  ! Arbitrary gamma is set using diagonal term of A
  u(l,l) = g(l)
  v(l,l) = 1.0_8 
end do

do j = 1, k
  a(j,j) = a(j,j) - g(j) 
  do i = n-k+j-1,n
    u(i,j) = m(i,j)
    v(i,j) = m(j,i)/g(j)
  end do
end do

do j = n-k+1,n
  do i = n-k+1,n
    a(i,j) = a(i,j) - sum(u(i,:)*v(j,:)) !This sum can be optimized
  end do
end do

ldab=2*k+k+1
allocate(ab(ldab,n))
ab = 0.0
do j = 1, n
  do i = max(1,j-k), min(n,j+k)
    ab(k+k+1+i-j,j) = a(i,j) 
  end do
end do
call print_matrix(ab)
call dgbtrf(n,n,k,k,ab,ldab,ipiv,info)

print*, "A:"; call print_matrix(a)
print*, 'U:'; call print_matrix(u)
print*, 'V:'; call print_matrix(v)
print*, 'M - (A + U.t(V)) : must be zero '
write(*,"(' Decomposition error = ', g15.3)") &
  sum(abs(m-(a+matmul(u,transpose(v)))))
call print_matrix(matmul(u,transpose(v)))

y=b
nrhs = 1
call dgbtrs('N',n,k,k,nrhs,ab,ldab,ipiv,y,n,info)
print*, ' Solve A.y = b error : ', sum(abs(b-matmul(a,y)))
nrhs = k
z = u
call dgbtrs('N',n,k,k,nrhs,ab,ldab,ipiv,z,n,info)
print*, ' Solve A.z = u error : ', sum(abs(u-matmul(a,z)))

!compute the matrix H = inverse(1+t(v).z)
call print_matrix(1.0+matmul(z,transpose(v)))
h = one + matmul(transpose(v),z)

print*,'H:';call print_matrix(h)
call dgetrf(k,k,h,k,jpiv,info)
call dgetri(k,h,k,jpiv,work,k*k,info)

print*, 'inverse(H)'
print*,'H:';call print_matrix(h)

print*, ' X = Y - Z . [H . (t(V).Y)] : '
x = y - matmul(z,matmul(h,matmul(transpose(v),y)))
call print_vector(x)
write(*,"(' error = ', g15.3)") sum(b - matmul(m,x))

!!store banded matrix A in q for banfac
!q = 0.0_8
!do j = 1, n
!  l = 0
!  do i = -k,k
!    l = l+1
!    if (i+j >= 1 .and. i+j <= n) q(l,j) = a(i+j,j)
!  end do
!end do
!
!write(*,*) "Banded matrix A stored in Q:"
!call print_matrix(q)
!
!!Factorize the matrix A
!call banfac ( q, k+kp1, n, k, k, iflag )
!print*, 'iflag=', iflag
!
!!Solve A.y = b
!x = b
!call banslv ( q, k+kp1, n, k, k, x )
!print*, ' Solve A.y = b error : ', sum(abs(b-matmul(a,x)))
!
!!Solve A.z = u
!do l = 1, k
!  call banslv ( q, k+kp1, n, k, k, u(:,l) )
!end do
!print*,'Z:';call print_matrix(u)
!!compute the matrix H = inverse(1+t(v).z)
!h = one+matmul(transpose(v),u)
!deallocate(ipiv); allocate(ipiv(k))
!
!call dgetrf(k,k,h,k,ipiv,info)
!call dgetri(k,h,k,ipiv,work,k*k,info)
!print*,'H:';call print_matrix(h)
!
!f = matmul(h,matmul(transpose(v),x))
!print*, ' X = Y - Z . [H . (t(V).Y)] : '
!x = x - matmul(u,f)
!call print_vector(x)
!write(*,"(' error = ', g15.3)") sum(b - matmul(m,x))

contains

subroutine print_vector(v)
real(8), intent(in)  :: v(:)
integer :: j
character(len=20) :: display_format

write(display_format, "('(''|''',i2,a,''' |'')')")size(v),'f9.3'
write(*,display_format) v
write(*,*)
end subroutine print_vector

subroutine print_matrix(m)
real(8), intent(in)  :: m(:,:)
integer :: j
character(len=20) :: display_format

write(display_format, "('(''|''',i2,a,''' |'')')")size(m,2),'f9.3'
do i = 1, size(m,1)
  write(*,display_format) m(i,:)
end do
write(*,*)
end subroutine print_matrix

subroutine outer_product(u,v,w)
real(8), intent(in)   :: u(:)
real(8), intent(in)   :: v(:)
real(8), intent(out)  :: w(:,:)

do j = 1, size(v)
  do i = 1, size(u)
    w(i,j) = u(i)*v(j)
  end do
end do
end subroutine outer_product

end program test_cyclic_banded_solver

