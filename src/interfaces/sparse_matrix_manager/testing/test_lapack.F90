program test_lapack
#include "sll_working_precision.h"
#include "sll_memory.h"

implicit none

sll_int32,  parameter   :: n = 5
sll_real64              :: A(n,n)
sll_real64              :: LU(n,n)
sll_real64              :: B(n)
sll_real64              :: X(n)
sll_real64, allocatable :: AB(:,:)
sll_int32,  allocatable :: ipiv(:)

sll_int32 :: kd, kl, ku, ldab, info
sll_int32 :: istep, nstep = 1
sll_int32 :: i, j, nrhs

A = 0.0_f64
A(1,1) = 2.0_f64; A(1,2) = -1.0_f64
do i = 2, n-1
   A(i,i  ) =  2.0_f64
   A(i,i-1) = -1.0_f64
   A(i,i+1) = -1.0_f64
end do
A(n,n)   = 2.0_f64
A(n,n-1) = -1.0_f64

do i = 1, n
   B(i) = i-1.0_f64
end do

do i = 1, n
   write(*,100)(sngl(A(i,j)), j= 1, n), B(i)
end do
write(*,*)

!General Matrix Factorization
allocate(IPIV(n))

LU = A
call DGETRF(n,n,LU,n,IPIV,INFO)
X = B
NRHS = 1
do istep=1, nstep
   call DGETRS('N',n,NRHS,LU,n,IPIV,X,n,INFO)
end do
write(*,100) X

!LU factorization of a real m-by-n band matrix
KL = 1
KU = 1
LDAB=2*KL+KU+1
allocate(AB(ldab,n))
AB = 0.0_f64
do j = 1, n
do i = max(1,j-KU), min(n,j+KL)
   AB(KL+KU+1+i-j,j) = A(i,j) 
end do
end do
call DGBTRF(n,n,KL,KU,AB,LDAB,IPIV,INFO)
X=B
NRHS = 1
do istep = 1, nstep
call DGBTRS('N',n,KL,KU,NRHS,AB,LDAB,IPIV,X,n,INFO)
end do
write(*,100) X

!Cholesky factorization of a real symmetric positive definite band matrix
KD = 1
LDAB=KD+1
deallocate(AB)
allocate(AB(ldab,n))
AB = 0.0_f64
do j = 1, n
do i = max(1,j-KD), j
   AB(KD+1+i-j,j) = A(i,j) 
end do
end do
NRHS = 1
call DPBTRF('U',n,KD,AB,LDAB,INFO)
X=B
do istep=1, nstep
   call DPBTRS('U',n,KD,NRHS,AB,LDAB,X,n,INFO)
end do
write(*,100) X

100 format(6(1x,f7.4))

print *,'#PASSED'

end program test_lapack
