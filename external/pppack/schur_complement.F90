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

  real(8), allocatable :: bb(:,:)
  real(8), allocatable :: cc(:,:)
  real(8), allocatable :: dd(:,:)
  real(8), allocatable :: yy(:,:)
  
  real(8), allocatable :: c1(:)
  real(8), allocatable :: c2(:)
  real(8), allocatable :: z1(:)
  real(8), allocatable :: z2(:)

end type schur_complement_solver

contains

subroutine schur_complement_fac(s, n, k, q)

  type(schur_complement_solver) :: s
  
  integer, intent(in)  :: n
  integer, intent(in)  :: k
  real(8)              :: q(2*k+1,n)
  integer              :: i
  integer              :: j
  integer              :: l
  integer              :: info
  integer              :: jpiv(k)
  integer              :: work(k*k)
  
  integer :: kp1
  
  kp1 = k+1
  allocate(s%bb(n-k,k  )); s%bb = 0.0_8
  allocate(s%cc(k  ,n-k)); s%cc = 0.0_8
  allocate(s%dd(k  ,k  )); s%dd = 0.0_8
  allocate(s%z2(n-k))    ; s%z2 = 0.0_8
  allocate(s%c2(k))      ; s%c2 = 0.0_8
  
  do i = 1, k
    l = 0
    do j = i, k
      s%bb(i,j) = q(k+kp1-l,n-k+l+i)
      l =l+1
      s%bb(n-k-k+j,l) = q(i,n-k+l)
    end do
  end do
  
  do j = 1, k
    l = 0
    do i = j, k
      l =l+1
      s%cc(i,j) = q(l,j)
    end do
  end do
  do i = 1, k
    l = 0
    do j = n-k-k+i,n-k
      s%cc(i,j) = q(kp1+k-l,n-k-k+l+i)
      l=l+1
    end do
  end do
  
  do i = 1, k
    l = 0
    do j = 1, k
      l=l+1
      s%dd(i,j) = q(kp1-l+i,n-k+l)
    end do
  end do
  
  write(*,*) "B"; call print_matrix(s%bb)
  write(*,*) "C"; call print_matrix(s%cc)
  write(*,*) "D"; call print_matrix(s%dd)
  
  !Factorize the matrix A
  call banfac ( q(:,1:n-k), k+kp1, n-k, k, k, info )
  
  !Solve A.Y = B
  allocate(s%yy(n-k,k))
  s%yy = s%bb
  do j = 1, k
    call banslv ( q(:,1:n-k), k+kp1, n-k, k, k, s%yy(:,j) )
  end do
  call print_matrix(s%yy)
  
  !Compute H= D - C.Y
  s%dd = s%dd - matmul(s%cc,s%yy)
  call print_matrix(s%dd)
  call dgetrf(k,k,s%dd,k,jpiv,info)
  call dgetri(k,s%dd,k,jpiv,work,k*k,info)
  

end subroutine schur_complement_fac

subroutine schur_complement_slv(s, n, k, q, x)

  type(schur_complement_solver) :: s
  integer, intent(in) :: n
  integer, intent(in) :: k
  real(8)             :: q(2*k+1,n)
  integer             :: kp1
  real(8)             :: x(n)
  
  kp1 = k+1
  
  !Solve A.z2 = b1
  s%z2 = x(1:n-k)
  call banslv ( q(:,1:n-k), k+kp1, n-k, k, k, s%z2 )
  !compute c2 = b2 - C.z2
  s%c2 = x(n-k+1:n) - matmul(s%cc,s%z2)
  !Solve H.x2 = c2
  x(n-k+1:n) = matmul(s%dd,s%c2)
  !Solve A.x1 = b1 - B.x2
  x(1:n-k) = x(1:n-k) - matmul(s%bb,x(n-k+1:n))
  call banslv ( q(:,1:n-k), k+kp1, n-k, k, k, x(1:n-k) )

end subroutine schur_complement_slv

subroutine print_vector(v)

  real(8), intent(in) :: v(:)
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


