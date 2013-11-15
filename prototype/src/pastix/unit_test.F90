program test_pastix
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_collective
use sll_pastix 
use sll_murge

implicit none

sll_int32                             :: npts = 100
sll_real64, dimension(:), allocatable :: sol
sll_real64, dimension(:), allocatable :: rhs
sll_int32                             :: error

call sll_boot_collective()

SLL_CLEAR_ALLOCATE(sol(1:npts), error)
SLL_CLEAR_ALLOCATE(rhs(1:npts), error)

call test_pastix_fortran()

!call test_pastix_murge()

call sll_halt_collective()

contains

subroutine test_pastix_fortran()
type(pastix_solver) :: solver

   call initialize_pastix(solver, npts)
   call solve_pastix(solver, sol, rhs)
   call delete_pastix(solver)

end subroutine test_pastix_fortran

subroutine test_pastix_murge()

   call initialize_murge()
   call delete_murge()

end subroutine test_pastix_murge

end program test_pastix
