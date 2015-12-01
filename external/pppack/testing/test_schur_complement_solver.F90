! To solve M.x = b


!Notre système peut s’écrire sous forme bloc
!(A B)
!(C D)
!où C et B sont des petites matrices. Le système devient
!A x1 + B x2 = b1
!C x1 + D x2 = b2
!On peut alors éliminer x1 de la deuxième equation en soustrayant CA^{-1} 
!fois la premiere:
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
!Resoudre: A Y = B  (Y=A^{-1}B), Y est une petite matrice. On peut appeler 
!dgbtrs  avec plusieurs second membres
!Calculer H = D - CY
!Factoriser H
!
!Solve:
!Calculer  c2=b2 -  C A^{-1} b1 (en résolvant  A z2=b1, puis b2=C z2) 
!Resoudre H x_2 = c2 puis A x1 = b1- B x2


program test_schur_complement_solver

use schur_complement

implicit none

integer, parameter :: n = 9
integer, parameter :: k = 3
integer, parameter :: kp1 = k+1

real(8) :: x(n)
real(8) :: b(n)
real(8) :: q(2*k+1,n)

real(8), allocatable :: aa(:,:)
real(8), allocatable :: bb(:,:)
real(8), allocatable :: cc(:,:)
real(8), allocatable :: dd(:,:)
real(8), allocatable :: yy(:,:)
real(8), allocatable :: hh(:,:)

real(8), allocatable :: x1(:)
real(8), allocatable :: x2(:)
real(8), allocatable :: b1(:)
real(8), allocatable :: b2(:)
real(8), allocatable :: c1(:)
real(8), allocatable :: c2(:)
real(8), allocatable :: z1(:)
real(8), allocatable :: z2(:)

integer :: i
integer :: j
integer :: l

real(8) :: m(n,n)
real(8) :: p(n,n)

integer              :: info
integer              :: nrhs
integer, allocatable :: ipiv(:)
integer              :: jpiv(k)
integer              :: work(k*k)

!Set the RHS 
do i = 1, n
  b(i) = real(i,8)
end do
!Create a cyclic banded system m
m = 0.0_8
do i = 1, n
  do j = -k,k
    l=modulo(j+i-1,n)+1 
    m(i,l) = real(i*10+l,8)
    m(l,i) = -real(l*10+i,8)
  end do
end do
write(*,*) "M="
call print_matrix(m)
!General Matrix Factorization with Lapack
allocate(ipiv(n)); ipiv=0
p = m
call dgetrf(n,n,p,n,ipiv,info)
x = b
nrhs = 1
call dgetrs('N',n,nrhs,p,n,ipiv,x,n,info)
write(*,"(' Lapack error = ', g15.3)") sum(abs(b-matmul(m,x)))
call print_vector(x)

q = 0.0_8
do j = 1, n
  l = 0
  do i = -k,k
    l = l+1
    q(l,j) = m(modulo(i+j-1,n)+1,j)
  end do
end do
write(*,*) "Q ="
call print_matrix(q)
write(*,*) "Banded matrix M stored in Q:"

x = b

call schur_complement_fac(n, k, q, x)

print*, ' x = '
call print_vector(x)
write(*,"(' error = ', g15.3)") sum(b - matmul(m,x))


end program test_schur_complement_solver


