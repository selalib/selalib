  
module sll_pastix

#define MPI_MASTER 0
#include "pastix_fortran.h"
#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_collective

implicit none
private

pastix_data_ptr_t                     :: pastix_data ! PaStiX structure (0 for first call)
pastix_int_t                          :: n           ! Number of columns in the matrix
pastix_int_t   ,dimension(:), pointer :: ia          ! Index of first elt of each col in ja and avals
pastix_int_t   ,dimension(:), pointer :: ja          ! Row of each element
pastix_float_t ,dimension(:), pointer :: avals       ! Value of each element
pastix_int_t   ,dimension(:), pointer :: perm        ! permutation tabular
pastix_int_t   ,dimension(:), pointer :: invp        ! reverse permutation tabular
pastix_float_t ,dimension(:), pointer :: rhs         ! Right hand side
pastix_int_t                          :: nrhs        ! right hand side number (only one possible)
pastix_int_t                          :: iparm(IPARM_SIZE) ! sll_int32 parameters
double precision                      :: dparm(DPARM_SIZE) ! Floating poin parameters
sll_int32                             :: nbthread    ! Number of threads in PaStiX
sll_int32                             :: verbose     ! Verbose level

sll_int32, private                    :: comm    
sll_int32, private                    :: prank    
sll_int32, private                    :: psize

public :: initialize_pastix, factorize_pastix

contains

subroutine initialize_pastix(npts)
sll_int32 :: npts
pastix_int_t          :: nnzero
sll_int32 :: i
sll_int32 :: j
sll_int32 :: error

n = npts
prank = sll_get_collective_rank( sll_world_collective )
psize = sll_get_collective_size( sll_world_collective )
comm  = sll_world_collective%comm

! Get options ftom command line

nbthread = 1
verbose  = API_VERBOSE_YES

nnzero =  3*n - 2

! Allocating
SLL_ALLOCATE( ia   (n+1)   , error)
SLL_ALLOCATE( ja   (nnzero), error)
SLL_ALLOCATE( avals(nnzero), error)
SLL_ALLOCATE( rhs  (n)     , error)

! Building ia, ja and avals and rhs
j=1
do i = 1, n
   ia(i) = j
   ! /* ONLY triangular inferior matrix */
   ja(j)    = i
   avals(j) = 2
   j=j+1
   if (i /= n) then
      ja(j)    = i+1
      avals(j) = -1.
      j=j + 1
   end if
   rhs(i) = 0
end do
ia(n+1) = j
rhs(1)  = 1
rhs(n)  = 1

! First PaStiX call to initiate parameters

pastix_data                   = 0
nrhs                          = 1
iparm(IPARM_MODIFY_PARAMETER) = API_NO
iparm(IPARM_START_TASK)       = API_TASK_INIT
iparm(IPARM_END_TASK)         = API_TASK_INIT

call pastix_fortran(pastix_data ,comm, &
                    n,ia,ja,avals,perm,invp,  &
                    rhs,nrhs,iparm,dparm)

! Customize some parameters

iparm(IPARM_THREAD_NBR) = nbthread  
iparm(IPARM_VERBOSE)    = verbose

iparm(IPARM_SYM)           = API_SYM_YES   ! API_SYM_NO
iparm(IPARM_FACTORIZATION) = API_FACT_LDLT ! API_FACT_LU

iparm(IPARM_MATRIX_VERIFICATION) = API_YES
iparm(IPARM_RHS_MAKING)          = API_RHS_1

SLL_ALLOCATE(perm(n), error)
SLL_ALLOCATE(invp(n), error)

! Call PaStiX first steps (Scotch - Fax - Blend)
iparm(IPARM_START_TASK) = API_TASK_ORDERING
iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

call pastix_fortran(pastix_data ,comm, &
                    n,ia,ja, &
                    avals,perm,invp,  &
                    rhs,nrhs,iparm,dparm)

end subroutine initialize_pastix

subroutine factorize_pastix()
    
! Call PaStiX factorization
iparm(IPARM_START_TASK) = API_TASK_NUMFACT
iparm(IPARM_END_TASK)   = API_TASK_NUMFACT
call pastix_fortran(pastix_data ,comm, &
                    n,ia,ja, &
                    avals,perm,invp,  &
                    rhs,nrhs,iparm,dparm)

comm = sll_world_collective%comm

! Call PaStiX updown and refinement
iparm(IPARM_START_TASK)       = API_TASK_SOLVE
iparm(IPARM_END_TASK)         = API_TASK_REFINE
! rhs will be changed to solution 
call pastix_fortran(pastix_data ,comm, &
                    n,ia,ja, &
                    avals,perm,invp,  &
                    rhs,nrhs,iparm,dparm)

comm = sll_world_collective%comm

! Call PaStiX clean
iparm(IPARM_START_TASK)       = API_TASK_CLEAN
iparm(IPARM_END_TASK)         = API_TASK_CLEAN
    
call pastix_fortran(pastix_data ,comm,         &
                    n,ia,ja,         &
                    avals,perm,invp, &
                    rhs,nrhs,iparm,dparm)
    
deallocate(ia)
deallocate(ja)
deallocate(avals)
deallocate(perm)
deallocate(invp)
deallocate(rhs)

end subroutine factorize_pastix
end module sll_pastix
