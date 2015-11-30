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


module schur_complement

implicit none

type :: schur_complement_solver

end type schur_complement_solver

contains

subroutine schur_complement_fac(n, k, q, m)

integer, intent(in) :: n
integer, intent(in) :: k
real(8)             :: q(2*k+1,n-k)
real(8)             :: m(n,n)

real(8) :: x(n)
real(8) :: b(n)

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

real(8) :: p(n,n)

integer              :: info
integer              :: nrhs
integer, allocatable :: ipiv(:)
integer              :: jpiv(k)
integer              :: work(k*k)

integer :: kp1

kp1 = k+1

do i = 1, n
  b(i) = real(i,8)
end do
allocate(b1(n-k)); b1 = b(1:n-k)
allocate(b2(k));   b2 = b(n-k+1:n)
allocate(x1(n-k)); x1 = 0.0_8
allocate(x2(k));   x2 = 0.0_8

!General Matrix Factorization with Lapack
allocate(ipiv(n)); ipiv=0
p = m
call dgetrf(n,n,p,n,ipiv,info)
x = b
nrhs = 1
call dgetrs('N',n,nrhs,p,n,ipiv,x,n,info)
write(*,"(' Lapack error = ', g15.3)") sum(abs(b-matmul(m,x)))
call print_vector(x)

allocate(aa(n-k,n-k))
aa = 0.0_8
do j = 1,n-k
  do i = max(1,j-k), min(n-k,j+k)
    aa(i,j) = m(i,j)
  end do
end do
write(*,*) "A="
call print_matrix(aa)

allocate(bb(n-k,k))
bb = 0.0_8
do j = 1, k
  do i = 1, n-k
    bb(i,j) = m(i,j+n-k)
  end do
end do
write(*,*) "B="
call print_matrix(bb)

allocate(cc(k,n-k))
cc = 0.0_8
do j = 1, n-k
  do i = 1, k
    cc(i,j) = m(i+n-k,j)
  end do
end do
write(*,*) "C="
call print_matrix(cc)

allocate(dd(k,k))
dd = 0.0_8
do j = 1, k
  do i = 1, k
    dd(i,j) = m(i+n-k,j+n-k)
  end do
end do
write(*,*) "D="
call print_matrix(dd)
!store banded matrix A in Q for banfac
q = 0.0_8
do j = 1, n-k
  l = 0
  do i = -k,k
    l = l+1
    if (i+j >=1 .and. i+j <=n-k ) q(l,j) = aa(i+j,j)
  end do
end do
write(*,*) "Q ="
call print_matrix(q)

write(*,*) "Banded matrix M stored in Q:"

!Factorize the matrix A
call banfac ( q, k+kp1, n-k, k, k, info )

!Solve A.Y = B
allocate(yy(n-k,k))
yy = bb
do j = 1, k
  call banslv ( q, k+kp1, n-k, k, k, yy(:,j) )
end do
write(*,*) "Y ="
call print_matrix(yy)

!Compute H= D - C.Y
allocate(hh(k,k))
hh = dd - matmul(cc,yy)
write(*,*) "H ="
call print_matrix(hh)
call dgetrf(k,k,hh,k,jpiv,info)
call dgetri(k,hh,k,jpiv,work,k*k,info)
print*,'H^(-1) =';call print_matrix(hh)

!Solve A.z2 = b1
allocate(z2(n-k))
z2 = b1
call banslv ( q, k+kp1, n-k, k, k, z2 )
write(*,*) " z2 = "
call print_vector(z2)

!compute c2 = b2 - C.z2
allocate(c2(k))
c2 = b2 - matmul(cc,z2)
write(*,*) " c2 = "
call print_vector(c2)
!Solve H.x2 = c2
x2 = matmul(hh,c2)
write(*,*) " x2 = "
call print_vector(x2)
!Solve A.x1 = b1 - B.x2

x1 = b1 - matmul(bb,x2)
call banslv ( q, k+kp1, n-k, k, k, x1 )
write(*,*) " x1 = "
call print_vector(x1)

x(1:n-k)   = x1
x(n-k+1:n) = x2

print*, ' x = '
call print_vector(x)
write(*,"(' error = ', g15.3)") sum(b - matmul(m,x))

end subroutine schur_complement_fac

subroutine print_vector(v)

  real(8), intent(in) :: v(:)
  integer             :: j
  character(len=20)   :: display_format
  
  write(display_format, "('(''|''',i2,a,''' |'')')")size(v),'f9.3'
  write(*,display_format) v
  write(*,*)

end subroutine print_vector

subroutine print_matrix(m)

  real(8), intent(in) :: m(:,:)
  integer             :: i
  character(len=20)   :: display_format
  
  write(display_format, "('(''|''',i2,a,''' |'')')")size(m,2),'f9.3'
  do i = 1, size(m,1)
    write(*,display_format) m(i,:)
  end do
  write(*,*)

end subroutine print_matrix

end module schur_complement


