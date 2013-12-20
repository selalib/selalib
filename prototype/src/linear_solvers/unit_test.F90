program test_linear_solvers
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_linear_solvers.h"

use sll_collective

implicit none

sll_real64, dimension(:), allocatable :: sol
sll_real64, dimension(:), allocatable :: sol_exact
sll_real64, dimension(:), allocatable :: rhs
type(sll_csc_matrix)                  :: csc_matrix

#ifdef PASTIX
type(pastix_solver)                   :: pastix
#endif

sll_int32   :: n, nnzeros
sll_real64  :: delta
sll_real64  :: x
sll_int32   :: error
sll_int32   :: i, j

call sll_boot_collective()

n = 32

SLL_ALLOCATE(sol(n),error)
SLL_ALLOCATE(rhs(n),error)
SLL_ALLOCATE(sol_exact(n),error)


delta = 1.0_f64 / n
do i=1,n
   x = (i-1)*delta
   rhs(i)       = 4*sll_pi**2*sin(2*sll_pi*x) * delta**2
   sol_exact(i) = sin(2*sll_pi*x)
end do
sol = rhs

nnzeros = 3*(n-1)+2

call initialize(csc_matrix, n-1, nnzeros, error)

j=1
do i = 1, n-1
   csc_matrix%colptr(i) = j
   csc_matrix%row(j)    = i
   csc_matrix%avals(j)  = 2
   j=j+1
   if (i /= n-1) then
      csc_matrix%row(j)   = i+1
      csc_matrix%avals(j) = -1.
      j=j+1
   end if
end do
csc_matrix%colptr(n) = j

!call initialize(pastix, n-1, nnzeros)
!pastix%colptr = csc_matrix%colptr
!pastix%row    = csc_matrix%row
!pastix%avals  = csc_matrix%avals
!call factorize(pastix)
!print*, sol
!print*, rhs
!print*, size(pastix%rhs), size(rhs), size(sol)
!stop
!call solve(pastix, sol(:))
!stop
!
!do i = 1, n
!   write(18,*) sngl((i-1)*delta), sngl(sol(i)), sngl(sol_exact(i))
!end do
!    
!print*,'   error=',maxval(abs(sol-sol_exact))

call sll_halt_collective()

end program test_linear_solvers
