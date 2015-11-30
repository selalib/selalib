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


program test_cyclic_banded_solver
implicit none

integer, parameter :: n = 9
integer, parameter :: k = 3

real(8) :: x(n)
real(8) :: y(n-k)
real(8) :: b(n)
real(8) :: q(2*k+1,n-k)
real(8) :: aa(n-k,n-k)
real(8) :: bb(k,k)
real(8) :: cc(k,k)
real(8) :: dd(k,k)

real(8) :: y(n-k)

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

integer       :: iflag

kp1 = k+1
do i = 1, n
  b(i) = 1.0_8 *i
end do

!Create a cyclic banded system m
m = 0.0_8
do i = 1, n
  do j = -k,k
    l=modulo(j+i-1,n)+1 
    m(i,l) = real(i*10+l,8)
  end do
end do
write(*,*) "M="
call print_matrix(m)

aa = 0.0_8
do j = 1,n-k
  do i = max(1,j-k), min(n-k,j+k)
    aa(i,j) = m(i,j)
  end do
end do
write(*,*) "A="
call print_matrix(aa)

bb = 0.0_8
do j = 1, k
  do i = 1, k
    bb(i,j) = m(i,j+n-k)
  end do
end do
write(*,*) "B="
call print_matrix(bb)

cc = 0.0_8
do j = 1, k
  do i = 1, k
    cc(i,j) = m(i+n-k,j)
  end do
end do
write(*,*) "C="
call print_matrix(cc)

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
call print_matrix(q)

write(*,*) "Banded matrix M stored in Q:"

!Factorize the matrix A
call banfac ( q, k+kp1, n-k, k, k, iflag )
print*, 'iflag=', iflag

!Solve A.Y = B
y = 0.0
do i = 1, k
  y = bb(1:n-k)
  call banslv ( q, k+kp1, n-k, k, k, x(1:n-k) )
end do
print*, ' Solve A.y = b error : ', sum(abs(b(1:n-k)-matmul(aa,x(1:n-k))))

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

  real(8), intent(in) :: v(:)
  integer             :: j
  character(len=20)   :: display_format
  
  write(display_format, "('(''|''',i2,a,''' |'')')")size(v),'f9.3'
  write(*,display_format) v
  write(*,*)

end subroutine print_vector

subroutine print_matrix(m)

  real(8), intent(in) :: m(:,:)
  integer             :: j
  character(len=20)   :: display_format
  
  write(display_format, "('(''|''',i2,a,''' |'')')")size(m,2),'f9.3'
  do i = 1, size(m,1)
    write(*,display_format) m(i,:)
  end do
  write(*,*)

end subroutine print_matrix

end program test_cyclic_banded_solver


