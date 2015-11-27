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

!store banded matrix A in q for banfac
q = 0.0_8
do j = 1, n
  l = 0
  do i = -k,k
    l = l+1
    q(l,j) = m(modulo(i+j-1,n)+1,j)
  end do
end do
call print_matrix(q)

call woodbury_solver(m, q, n, k, b)

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


