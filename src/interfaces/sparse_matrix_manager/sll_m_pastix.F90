module sll_m_pastix
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define MPI_MASTER 0
#include "pastix_fortran.h"
#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_m_collective, only:  &
  sll_s_boot_collective,     &
  sll_f_get_collective_rank, &
  sll_f_get_collective_size, &
  sll_v_world_collective

implicit none

public :: &
  pastix_solver, &
  initialize, &
  factorize, &
  solve, &
  delete

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
type :: pastix_solver
 pastix_data_ptr_t        :: pastix_data !< PaStiX structure (0 for first call)
 pastix_int_t             :: ncols       !< Number of columns in the matrix
 pastix_int_t   , pointer :: colptr(:)   !< Index of first elt of each col in avals
 pastix_int_t   , pointer :: row(:)      !< Row of each element
 pastix_float_t , pointer :: avals(:)    !< Value of each element
 pastix_int_t   , pointer :: perm(:)     !< permutation tabular
 pastix_int_t   , pointer :: invp(:)     !< reverse permutation tabular
 pastix_float_t , pointer :: rhs(:)      !< Right hand side
 pastix_int_t             :: nrhs        !< right hand side number (only one possible)
 pastix_int_t             :: iparm(IPARM_SIZE) ! sll_int32 parameters
 double precision         :: dparm(DPARM_SIZE) ! Floating poin parameters
 sll_int32                :: nbthread    !< Number of threads in PaStiX
 sll_int32                :: verbose     !< Verbose level
end type pastix_solver


interface initialize
   module procedure init_pastix
end interface initialize

interface solve
   module procedure solve_pastix_with_rhs
   module procedure solve_pastix_without_rhs
end interface solve

interface factorize
   module procedure factorize_pastix
end interface factorize

interface delete
   module procedure free_pastix
end interface delete

contains

subroutine init_pastix(self,n,nnzeros,row_ptr,col_ind,val)

  type(pastix_solver)              :: self
  sll_int32,    intent(in)         :: n
  pastix_int_t, intent(in)         :: nnzeros
  sll_int32,    intent(in), target :: row_ptr(:)
  sll_int32,    intent(in), target :: col_ind(:)
  sll_real64,   intent(in), target :: val(:)
  sll_int32                        :: error
  sll_int32                        :: comm    
  sll_int32                        :: prank    
  sll_int32                        :: psize

  if( .not. associated(sll_v_world_collective)) then
     call sll_s_boot_collective()
  end if
  self%ncols = n
  prank = sll_f_get_collective_rank( sll_v_world_collective )
  psize = sll_f_get_collective_size( sll_v_world_collective )
  comm  = sll_v_world_collective%comm

  ! Get options ftom command line
  self%nbthread = 1
  self%verbose  = API_VERBOSE_NO

  ! Allocating
  self%colptr => row_ptr
  self%row    => col_ind
  self%avals  => val
  SLL_CLEAR_ALLOCATE(self%rhs(nnzeros), error)

  ! First PaStiX call to initiate parameters

  self%pastix_data                   = 0
  self%nrhs                          = 1
  self%iparm(IPARM_MODIFY_PARAMETER) = API_NO
  self%iparm(IPARM_START_TASK)       = API_TASK_INIT
  self%iparm(IPARM_END_TASK)         = API_TASK_INIT

  SLL_ALLOCATE(self%perm(n), error)
  SLL_ALLOCATE(self%invp(n), error)

  call pastix_fortran(self%pastix_data,comm,n,self%colptr,self%row, &
                      self%avals,self%perm,self%invp,self%rhs,    &
                      self%nrhs,self%iparm,self%dparm)

  ! Customize some parameters

  self%iparm(IPARM_THREAD_NBR) = self%nbthread  
  self%iparm(IPARM_VERBOSE)    = self%verbose

  self%iparm(IPARM_SYM)           = API_SYM_YES !API_SYM_NO
  self%iparm(IPARM_FACTORIZATION) = API_FACT_LU !API_FACT_LDLT

  self%iparm(IPARM_MATRIX_VERIFICATION) = API_YES
  
  !The matrix is in CSR format, we transpose to get CSC
  self%iparm(IPARM_TRANSPOSE_SOLVE) = API_YES
  
end subroutine init_pastix

subroutine factorize_pastix(self)

  type(pastix_solver) :: self
  sll_int32           :: comm    
   
  comm  = sll_v_world_collective%comm
  ! Call PaStiX first steps (Scotch - Fax - Blend)
  self%iparm(IPARM_START_TASK) = API_TASK_ORDERING
  self%iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

  call pastix_fortran(self%pastix_data, &
                      comm,             &
                      self%ncols,       &
                      self%colptr,      &
                      self%row,         &
                      self%avals,       &
                      self%perm,        &
                      self%invp,        &
                      self%rhs,         &
                      self%nrhs,        &
                      self%iparm,       &
                      self%dparm        )

  ! Call PaStiX factorization
  self%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
  self%iparm(IPARM_END_TASK)   = API_TASK_NUMFACT

  call pastix_fortran(self%pastix_data, &
                      comm,             &
                      self%ncols,       &
                      self%colptr,      &
                      self%row,         &
                      self%avals,       &
                      self%perm,        &
                      self%invp,        &
                      self%rhs,         &
                      self%nrhs,        &
                      self%iparm,       &
                      self%dparm)

end subroutine factorize_pastix

subroutine solve_pastix_without_rhs(self, sol)

  type(pastix_solver)              :: self
  sll_real64, dimension(:), target :: sol
  sll_int32                        :: comm    

  comm = sll_v_world_collective%comm

  self%rhs => sol

  ! Call PaStiX updown and refinement
  self%iparm(IPARM_START_TASK) = API_TASK_SOLVE
  self%iparm(IPARM_END_TASK)   = API_TASK_REFINE
  ! rhs will be changed to solution 
  call pastix_fortran(self%pastix_data ,comm,                  &
                      self%ncols,self%colptr,self%row,         &
                      self%avals,self%perm,self%invp,          &
                      self%rhs,self%nrhs,self%iparm,self%dparm)

  sol = self%rhs

end subroutine solve_pastix_without_rhs

subroutine solve_pastix_with_rhs(self, rhs, sol)

  type(pastix_solver)      :: self
  sll_real64, dimension(:) :: rhs
  sll_real64, dimension(:) :: sol
  sll_int32                :: comm    

  comm = sll_v_world_collective%comm

  self%rhs = rhs

  ! Call PaStiX updown and refinement
  self%iparm(IPARM_START_TASK) = API_TASK_SOLVE
  self%iparm(IPARM_END_TASK)   = API_TASK_REFINE
  ! rhs will be changed to solution 
  call pastix_fortran(self%pastix_data ,comm,                  &
                      self%ncols,self%colptr,self%row,         &
                      self%avals,self%perm,self%invp,          &
                      self%rhs,self%nrhs,self%iparm,self%dparm)

  sol = self%rhs

end subroutine solve_pastix_with_rhs

subroutine free_pastix(self)

  type(pastix_solver) :: self
  sll_int32           :: comm    

  comm = sll_v_world_collective%comm

  ! Call PaStiX clean
  self%iparm(IPARM_START_TASK)       = API_TASK_CLEAN
  self%iparm(IPARM_END_TASK)         = API_TASK_CLEAN
   
  call pastix_fortran(self%pastix_data ,comm,                  &
                      self%ncols,self%colptr,self%row,         &
                      self%avals,self%perm,self%invp,          &
                      self%rhs,self%nrhs,self%iparm,self%dparm)
   
  deallocate(self%colptr)
  deallocate(self%row)
  deallocate(self%avals)
  deallocate(self%perm)
  deallocate(self%invp)
  deallocate(self%rhs)
   
end subroutine free_pastix

end module sll_m_pastix
