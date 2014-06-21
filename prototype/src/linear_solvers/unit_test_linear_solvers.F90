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

sll_int32   :: ncells
sll_int32   :: n, nnzeros
sll_real64  :: delta
sll_real64  :: x
sll_int32   :: error
sll_int32   :: i
sll_int32   :: j

#ifdef PASTIX
call sll_boot_collective()
#endif

ncells = 33

SLL_CLEAR_ALLOCATE(sol(1:ncells),error)
SLL_CLEAR_ALLOCATE(rhs(1:ncells),error)
SLL_CLEAR_ALLOCATE(sol_exact(1:ncells),error)

delta = 1.0_f64 / ncells
do i=1,ncells
   x = (i-1)*delta
   rhs(i)       = 4*sll_pi**2*sin(2*sll_pi*x) * delta**2
   sol_exact(i) = sin(2*sll_pi*x)
end do
sol = rhs

n = ncells-1
nnzeros = 3*n+2

call initialize(csc_matrix, n, nnzeros, error)

#ifdef PASTIX
call initialize(pastix, n, nnzeros)

j=1
do i = 1, n
   pastix%colptr(i) = j
   pastix%row(j)    = i
   pastix%avals(j) = 2
   j=j+1
   if (i /= n) then
      pastix%row(j)   = i+1
      pastix%avals(j) = -1.
      j=j+1
   end if
end do
pastix%colptr(n+1) = j


!pastix%colptr = csc_matrix%colptr
!pastix%row    = csc_matrix%row
!pastix%avals  = csc_matrix%avals
call factorize(pastix)
!print*, sol
!print*, rhs
!print*, size(pastix%rhs), size(rhs), size(sol)
!stop
call solve(pastix, sol(2:ncells))
!stop

do i = 1, ncells
   write(18,*) sngl((i-1)*delta), sngl(sol(i)), sngl(sol_exact(i))
end do
    
print*,'   error=',maxval(abs(sol-sol_exact))


call sll_halt_collective()
#endif

end program test_linear_solvers
