  
module sll_pastix

#define MPI_MASTER 0
#include "pastix_fortran.h"
#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_collective

implicit none

type, public :: pastix_solver
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
end type pastix_solver

sll_int32, private                    :: comm    
sll_int32, private                    :: prank    
sll_int32, private                    :: psize

contains

subroutine initialize_pastix(this, n)
type(pastix_solver)   :: this
sll_int32, intent(in) :: n
pastix_int_t          :: nnzero
sll_int32 :: i
sll_int32 :: j
sll_int32 :: error

prank = sll_get_collective_rank( sll_world_collective )
psize = sll_get_collective_size( sll_world_collective )
comm  = sll_world_collective%comm

! Get options ftom command line

this%n = n
this%nbthread = 1
this%verbose  = API_VERBOSE_YES

nnzero =  3*n - 2

! Allocating
SLL_ALLOCATE( this%ia   (n+1)   , error)
SLL_ALLOCATE( this%ja   (nnzero), error)
SLL_ALLOCATE( this%avals(nnzero), error)
SLL_ALLOCATE( this%rhs  (n)     , error)
stop

! Building ia, ja and avals and rhs
j=1
do i = 1, n
   this%ia(i) = j
   ! /* ONLY triangular inferior matrix */
   this%ja(j)    = i
   this%avals(j) = 2
   j=j+1
   if (i /= this%n) then
      this%ja(j)    = i+1
      this%avals(j) = -1.
      j=j + 1
   end if
   this%rhs(i) = 0
end do
this%ia(n+1) = j
this%rhs(1)  = 1
this%rhs(n)  = 1

! First PaStiX call to initiate parameters

this%pastix_data                   = 0
this%nrhs                          = 1
this%iparm(IPARM_MODIFY_PARAMETER) = API_NO
this%iparm(IPARM_START_TASK)       = API_TASK_INIT
this%iparm(IPARM_END_TASK)         = API_TASK_INIT

call pastix_fortran(this%pastix_data ,comm, &
                    this%n,this%ia,this%ja,this%avals,this%perm,this%invp,  &
                    this%rhs,this%nrhs,this%iparm,this%dparm)

! Customize some parameters

this%iparm(IPARM_THREAD_NBR) = this%nbthread  
this%iparm(IPARM_VERBOSE)    = this%verbose

this%iparm(IPARM_SYM)           = API_SYM_YES   ! API_SYM_NO
this%iparm(IPARM_FACTORIZATION) = API_FACT_LDLT ! API_FACT_LU

this%iparm(IPARM_MATRIX_VERIFICATION) = API_YES
this%iparm(IPARM_RHS_MAKING)          = API_RHS_1

SLL_ALLOCATE(this%perm(n), error)
SLL_ALLOCATE(this%invp(n), error)

! Call PaStiX first steps (Scotch - Fax - Blend)
this%iparm(IPARM_START_TASK) = API_TASK_ORDERING
this%iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

call pastix_fortran(this%pastix_data ,comm, &
                    this%n,this%ia,this%ja, &
                    this%avals,this%perm,this%invp,  &
                    this%rhs,this%nrhs,this%iparm,this%dparm)
    
! Call PaStiX factorization
this%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
this%iparm(IPARM_END_TASK)   = API_TASK_NUMFACT
call pastix_fortran(this%pastix_data ,comm, &
                    this%n,this%ia,this%ja, &
                    this%avals,this%perm,this%invp,  &
                    this%rhs,this%nrhs,this%iparm,this%dparm)

end subroutine initialize_pastix

subroutine solve_pastix(this)
type(pastix_solver)      :: this

comm = sll_world_collective%comm

! Call PaStiX updown and refinement
this%iparm(IPARM_START_TASK)       = API_TASK_SOLVE
this%iparm(IPARM_END_TASK)         = API_TASK_REFINE
! rhs will be changed to solution 
call pastix_fortran(this%pastix_data ,comm, &
                    this%n,this%ia,this%ja, &
                    this%avals,this%perm,this%invp,  &
                    this%rhs,this%nrhs,this%iparm,this%dparm)

end subroutine solve_pastix

subroutine delete_pastix(this)
type(pastix_solver) :: this

comm = sll_world_collective%comm

! Call PaStiX clean
this%iparm(IPARM_START_TASK)       = API_TASK_CLEAN
this%iparm(IPARM_END_TASK)         = API_TASK_CLEAN
    
call pastix_fortran(this%pastix_data ,comm,         &
                    this%n,this%ia,this%ja,         &
                    this%avals,this%perm,this%invp, &
                    this%rhs,this%nrhs,this%iparm,this%dparm)
    
deallocate(this%ia)
deallocate(this%ja)
deallocate(this%avals)
deallocate(this%perm)
deallocate(this%invp)
deallocate(this%rhs)

end subroutine delete_pastix

end module sll_pastix
