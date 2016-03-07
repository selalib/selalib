#define MPI_MASTER 0

program test_pastix
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_m_collective
use sll_m_pastix 

implicit none

sll_int32 :: nhinv, m

type(pastix_solver)               :: linear_solver
sll_int32                         :: nnzeros
sll_int32,  dimension(:), pointer :: ia
sll_int32,  dimension(:), pointer :: ja
sll_real64, dimension(:), pointer :: avals
sll_real64, dimension(:), pointer :: rhs
sll_real64, dimension(:), pointer :: sol

call sll_s_boot_collective()

!set inverse of the mesh size (feel free to change)
nhinv=5
!first allocate the vectors with correct size
m=(nhinv-1)**2
nnzeros = 5*m
allocate(avals(nnzeros),ja(nnzeros),ia(m+1),rhs(m),sol(m))

call uni_laplace_2d(nhinv-1,rhs,avals,ja,ia)

call initialize(linear_solver,m,nnzeros,ia,ja,avals)

print *,'#enter factorize'
call factorize(linear_solver)
print *,'#end of factorize'
sol = rhs
print *,'#enter solve'
call solve(linear_solver,sol)
print *,'#end of solve'
write(*,100) sol
call delete(linear_solver)
print *,'#end of delete'

call sll_s_halt_collective()

print *,'#PASSED'
100 format(6(1x,f7.4))

contains

!----------------------------------------------------------------------
subroutine uni_laplace_2d(m,f,a,ja,ia)
!
! Fill a matrix in CSR format corresponding to a constant coefficient
! five-point stencil on a square grid
!
sll_real64 :: f(:)
sll_real64 :: a(:)
sll_int32  :: m
sll_int32  :: ia(:)
sll_int32  :: ja(:)
sll_int32  :: k,l,i,j

sll_real64, parameter :: zero =  0.0_f64
sll_real64, parameter :: cx   = -1.0_f64
sll_real64, parameter :: cy   = -1.0_f64
sll_real64, parameter :: cd   =  4.0_f64

k     = 0
l     = 0
ia(1) = 1
do i=1,m
  do j=1,m
    k     = k+1
    l     = l+1
    a(l)  = cd
    ja(l) = k
    f(k)  = zero
    if (j < m) then
      l     = l+1
      a(l)  = cx
      ja(l) = k+1
    else
      f(k) = f(k)-cx
    end if
    if (i < m) then
      l     = l+1
      a(l)  = cy
      ja(l) = k+m
    else
      f(k)  = f(k)-cy
    end if
    if (j > 1) then
      l     = l+1
      a(l)  = cx
      ja(l) = k-1
    else
      f(k)  = f(k)-cx
    end if
    if (i > 1) then
      l     = l+1
      a(l)  = cy
      ja(l) = k-m
    else
      f(k)  = f(k)-cy
    end if
    ia(k+1) = l+1
  end do
end do

end subroutine uni_laplace_2D

end program test_pastix
