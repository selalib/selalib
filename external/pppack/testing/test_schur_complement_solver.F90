program test_schur_complement_solver

use schur_complement

implicit none

type(schur_complement_solver) :: s
integer, parameter :: n = 9
integer, parameter :: k = 3

real(8) :: x(n)
real(8) :: b(n)
real(8) :: q(2*k+1,n)

integer :: i
integer :: j
integer :: l

real(8) :: m(n,n)

!Set the RHS 
do i = 1, n
  b(i) = real(i,8)
end do

!Create an arbitrary cyclic banded system m
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

q = 0.0_8
do j = 1, n
  l = 0
  do i = -k,k
    l = l+1
    q(l,j) = m(modulo(i+j-1,n)+1,j)
  end do
end do
write(*,*) "Banded matrix M stored in Q:"
call print_matrix(q)

x = b

call schur_complement_fac(s, n, k, q)
call schur_complement_slv(s, n, k, q, x)

print*, ' x = '; call print_vector(x)
write(*,"(' error = ', g15.3)") maxval(b - matmul(m,x))


end program test_schur_complement_solver


