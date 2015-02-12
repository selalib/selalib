program pastix_f90

use mpi
#define MPI_MASTER 0
#include "pastix_fortran.h"

implicit none

pastix_data_ptr_t                    :: pastix_data! PaStiX structure (0 for first call)
pastix_int_t                         :: n          ! Number of columns in the matrix
pastix_int_t  ,dimension(:), pointer :: ia         ! Index of first elt of each col in ja and avals
pastix_int_t  ,dimension(:), pointer :: ja         ! Row of each element
pastix_float_t,dimension(:), pointer :: avals      ! Value of each element
pastix_int_t  ,dimension(:), pointer :: perm       ! permutation tabular
pastix_int_t  ,dimension(:), pointer :: invp       ! reverse permutation tabular
pastix_float_t,dimension(:), pointer :: rhs        ! Right hand side
pastix_int_t                         :: nrhs       ! right hand side number (1)
pastix_int_t  ,dimension(IPARM_SIZE) :: iparm      ! integer parameters
real(8),       dimension(DPARM_SIZE) :: dparm      ! Floating poin parameters
pastix_int_t                         :: nnzero
integer                              :: nbthread   ! Number of threads in PaStiX
integer                              :: verbose    ! Verbose level
integer                              :: comm    
integer                              :: i
integer                              :: error
integer                              :: provided
integer                              :: required

!
! Initiate MPI communication
!
comm = MPI_COMM_WORLD
    
required=MPI_THREAD_MULTIPLE
call MPI_Init_thread(required,provided,error)

n = 5

! Get options ftom command line

nbthread = 1
verbose  = API_VERBOSE_NO
nnzero   =  3*n - 2

! Allocating
allocate( ia   (n+1)   ); ia = 0
allocate( ja   (nnzero)); ja = 0
allocate( avals(nnzero)); avals = 0.0
allocate( rhs  (n)     ); rhs = 0.0

ia(1:6) = (/1,3,5,7,9,10/)
ja(1:9) = (/1,2,2,3,3,4,4,5,5/)
avals(1:9) = (/2,-1,2,-1,2,-1,2,-1,2/)

write(*,"(a10,6i4)")   "ia    : ", ia(1:6)
write(*,"(a10,9i3)")   "ja    : ", ja(1:9)
write(*,"(a10,9f6.1)") "avals : ", avals(1:9)

! First PaStiX call to initiate parameters

pastix_data                   = 0
nrhs                          = 1
iparm(IPARM_MODIFY_PARAMETER) = API_NO
iparm(IPARM_START_TASK)       = API_TASK_INIT
iparm(IPARM_END_TASK)         = API_TASK_INIT

call pastix_fortran(pastix_data,comm,n,ia,ja, &
                    avals,perm,invp,rhs,nrhs,iparm,dparm)

! Customize some parameters

iparm(IPARM_THREAD_NBR)          = nbthread  
iparm(IPARM_VERBOSE)             = verbose

iparm(IPARM_SYM)                 = API_SYM_YES   ! API_SYM_NO
iparm(IPARM_FACTORIZATION)       = API_FACT_LDLT ! API_FACT_LU

iparm(IPARM_MATRIX_VERIFICATION) = API_YES

allocate(perm(n))
allocate(invp(n))

! Call PaStiX first steps (Scotch - Fax - Blend)
iparm(IPARM_START_TASK) = API_TASK_ORDERING
iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

call pastix_fortran(pastix_data,comm,n,ia,ja, &
                    avals,perm,invp,rhs,nrhs,iparm,dparm)

do i = 1, n
   rhs(i) = i-1
end do
! Call PaStiX factorization
iparm(IPARM_START_TASK) = API_TASK_NUMFACT
iparm(IPARM_END_TASK)   = API_TASK_NUMFACT
call pastix_fortran(pastix_data,comm,n,ia,ja, &
                    avals,perm,invp,rhs,nrhs,iparm,dparm)

do i = 1, n
   rhs(i) = i-1
end do
write(*,"(a,5f7.2)") " RHS : " , rhs

! Call PaStiX updown and refinement
iparm(IPARM_START_TASK) = API_TASK_SOLVE
iparm(IPARM_END_TASK)   = API_TASK_REFINE
! rhs will be changed to solution 
call pastix_fortran(pastix_data ,comm, &
                    n,ia,ja, &
                    avals,perm,invp,  &
                    rhs,nrhs,iparm,dparm)

write(*,"(a,5f7.2)") " SOL : " , rhs

! Call PaStiX clean
iparm(IPARM_START_TASK)       = API_TASK_CLEAN
iparm(IPARM_END_TASK)         = API_TASK_CLEAN
    
call pastix_fortran(pastix_data,comm,n,ia,ja, &
                    avals,perm,invp,rhs,nrhs,iparm,dparm)
    
deallocate(ia)
deallocate(ja)
deallocate(avals)
deallocate(perm)
deallocate(invp)
deallocate(rhs)

end program pastix_f90
