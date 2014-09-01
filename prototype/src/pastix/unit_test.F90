#define MPI_MASTER 0

program test_pastix
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_collective
use sll_pastix 
use sll_murge

implicit none

sll_int32, parameter                  :: NPTS = 5
sll_real64 :: A(NPTS,NPTS), LU(NPTS,NPTS), B(NPTS), X(NPTS)
sll_real64, allocatable :: AB(:,:)
sll_int32, allocatable :: ipiv(:)
sll_int32 :: kd, kl, ku, ldab, info
sll_int32 :: istep, nstep = 1
sll_int32 :: nrhs
type(pastix_solver) :: linear_solver

call sll_boot_collective()

A = 0
A(1,1) = 2; A(1,2) = -1.0
do i = 2, NPTS-1
   A(i,i) =  2
   A(i,i-1) = -1
   A(i,i+1) = -1
end do
A(NPTS,NPTS) = 2
A(NPTS,NPTS-1) = -1

do i = 1, NPTS
   B(i) = i-1
end do

do i = 1, NPTS
   write(*,100)(sngl(A(i,j)), j= 1, NPTS), B(i)
end do

!General Matrix Factorization
allocate(IPIV(NPTS))

LU = A
call DGETRF(NPTS,NPTS,LU,NPTS,IPIV,INFO)
X = B
NRHS = 1
do istep=1, nstep
   call DGETRS('N',NPTS,NRHS,LU,NPTS,IPIV,X,NPTS,INFO)
end do
write(*,100) X

!LU factorization of a real m-by-n band matrix
KL = 1
KU = 1
LDAB=2*KL+KU+1
allocate(AB(ldab,NPTS))
AB = 0.0
do j = 1, NPTS
do i = max(1,j-KU), min(NPTS,j+KL)
   AB(KL+KU+1+i-j,j) = A(i,j) 
end do
end do
call DGBTRF(NPTS,NPTS,KL,KU,AB,LDAB,IPIV,INFO)
X=B
NRHS = 1
do istep = 1, nstep
call DGBTRS('N',NPTS,KL,KU,NRHS,AB,LDAB,IPIV,X,NPTS,INFO)
end do
write(*,100) X


!Cholesky factorization of a real symmetric positive definite band matrix
KD = 1
LDAB=KD+1
deallocate(AB)
allocate(AB(ldab,NPTS))
AB = 0.0
do j = 1, NPTS
do i = max(1,j-KD), j
   AB(KD+1+i-j,j) = A(i,j) 
end do
end do
NRHS = 1
call DPBTRF('U',NPTS,KD,AB,LDAB,INFO)
X=B
do istep=1, nstep
   call DPBTRS('U',NPTS,KD,NRHS,AB,LDAB,X,NPTS,INFO)
end do
write(*,100) X

100 format(6(1x,f7.4))
!200 format(5(1x,f7.4))

call initialize(linear_solver,NPTS,3*NPTS-2)

linear_solver%colptr(1:6)    = (/1,3,5,7,9,10/)
linear_solver%row(1:9)    = (/1,2,2,3,3,4,4,5,5/)
linear_solver%avals(1:9) = (/2,-1,2,-1,2,-1,2,-1,2/)

write(*,"(a10,6i4)")   "ia    : ", linear_solver%colptr(1:6)
write(*,"(a10,9i3)")   "ja    : ", linear_solver%row(1:9)
write(*,"(a10,9f6.1)") "avals : ", linear_solver%avals(1:9)

print *,'#enter factorize'
call factorize(linear_solver)
print *,'#end of factorize'
X = B
print *,'#enter solve'
call solve(linear_solver,X)
print *,'#end of solve'
write(*,100) X
call delete(linear_solver)
print *,'#end of delete'


!call initialize_murge()
!call delete_murge()

call sll_halt_collective()

print *,'#PASSED'

end program test_pastix
