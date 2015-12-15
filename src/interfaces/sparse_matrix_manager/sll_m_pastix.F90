module sll_m_pastix
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define MPI_MASTER 0
#include "pastix_fortran.h"
#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_m_collective, only: &
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
   module procedure initialize_pastix
   !module procedure initialize_pastix_with_csc_matrix
end interface initialize

interface solve
   module procedure solve_pastix_with_rhs
   module procedure solve_pastix_without_rhs
end interface solve

interface factorize
   module procedure factorize_pastix
end interface factorize

interface delete
   module procedure delete_pastix
end interface delete

contains

subroutine initialize_pastix(this,n,nnzeros,row_ptr,col_ind,val)

  type(pastix_solver)                :: this
  sll_int32,    intent(in)           :: n
  pastix_int_t, intent(in)           :: nnzeros
  sll_int32,    intent(in), optional :: row_ptr(:)
  sll_int32,    intent(in), optional :: col_ind(:)
  sll_real64,   intent(in), optional :: val(:)
  sll_int32                          :: error
  sll_int32                          :: comm    
  sll_int32                          :: prank    
  sll_int32                          :: psize

  if( .not. associated(sll_v_world_collective)) then
     call sll_s_boot_collective()
  end if
  this%ncols = n
  prank = sll_f_get_collective_rank( sll_v_world_collective )
  psize = sll_f_get_collective_size( sll_v_world_collective )
  comm  = sll_v_world_collective%comm

  ! Get options ftom command line
  this%nbthread = 1
  this%verbose  = API_VERBOSE_NO

  ! Allocating
  SLL_ALLOCATE( this%colptr(1:n+1), error)
  if (present(row_ptr)) this%colptr = row_ptr
  SLL_ALLOCATE( this%row(1:nnzeros), error)
  if (present(col_ind)) this%row    = col_ind
  SLL_ALLOCATE( this%avals(1:nnzeros), error)
  if (present(val)) this%avals  = val
  SLL_ALLOCATE( this%rhs(1:n), error)

  ! First PaStiX call to initiate parameters

  this%pastix_data                   = 0
  this%nrhs                          = 1
  this%iparm(IPARM_MODIFY_PARAMETER) = API_NO
  this%iparm(IPARM_START_TASK)       = API_TASK_INIT
  this%iparm(IPARM_END_TASK)         = API_TASK_INIT

  SLL_ALLOCATE(this%perm(n), error)
  SLL_ALLOCATE(this%invp(n), error)

  call pastix_fortran(this%pastix_data,comm,n,this%colptr,this%row, &
                      this%avals,this%perm,this%invp,this%rhs,    &
                      this%nrhs,this%iparm,this%dparm)

  ! Customize some parameters

  this%iparm(IPARM_THREAD_NBR) = this%nbthread  
  this%iparm(IPARM_VERBOSE)    = this%verbose

  this%iparm(IPARM_SYM)           = API_SYM_YES !API_SYM_NO
  this%iparm(IPARM_FACTORIZATION) = API_FACT_LU !API_FACT_LDLT

  this%iparm(IPARM_MATRIX_VERIFICATION) = API_YES
  
  !The matrix is in CSR format, we transpose to get CSC
  this%iparm(IPARM_TRANSPOSE_SOLVE) = API_YES
  
end subroutine initialize_pastix

subroutine factorize_pastix(this)

  type(pastix_solver) :: this
  sll_int32           :: comm    
   
  comm  = sll_v_world_collective%comm
  ! Call PaStiX first steps (Scotch - Fax - Blend)
  this%iparm(IPARM_START_TASK) = API_TASK_ORDERING
  this%iparm(IPARM_END_TASK)   = API_TASK_ANALYSE

  call pastix_fortran(this%pastix_data,comm,this%ncols,this%colptr,     &
                      this%row,this%avals,this%perm,this%invp,    &
                      this%rhs,this%nrhs,this%iparm,this%dparm)

  ! Call PaStiX factorization
  this%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
  this%iparm(IPARM_END_TASK)   = API_TASK_NUMFACT
  call pastix_fortran(this%pastix_data,comm,this%ncols,this%colptr,this%row, &
                      this%avals,this%perm,this%invp,this%rhs,this%nrhs, &
                      this%iparm,this%dparm)

end subroutine factorize_pastix

subroutine solve_pastix_without_rhs(this, sol)

  type(pastix_solver)      :: this
  sll_real64, dimension(:) :: sol
  sll_int32                :: comm    

  comm = sll_v_world_collective%comm

  this%rhs = sol

  ! Call PaStiX updown and refinement
  this%iparm(IPARM_START_TASK) = API_TASK_SOLVE
  this%iparm(IPARM_END_TASK)   = API_TASK_REFINE
  ! rhs will be changed to solution 
  call pastix_fortran(this%pastix_data ,comm,                  &
                      this%ncols,this%colptr,this%row,         &
                      this%avals,this%perm,this%invp,          &
                      this%rhs,this%nrhs,this%iparm,this%dparm)

  sol = this%rhs

end subroutine solve_pastix_without_rhs

subroutine solve_pastix_with_rhs(this, rhs, sol)

  type(pastix_solver)      :: this
  sll_real64, dimension(:) :: rhs
  sll_real64, dimension(:) :: sol
  sll_int32                :: comm    

  comm = sll_v_world_collective%comm

  this%rhs = rhs

  ! Call PaStiX updown and refinement
  this%iparm(IPARM_START_TASK) = API_TASK_SOLVE
  this%iparm(IPARM_END_TASK)   = API_TASK_REFINE
  ! rhs will be changed to solution 
  call pastix_fortran(this%pastix_data ,comm,                  &
                      this%ncols,this%colptr,this%row,         &
                      this%avals,this%perm,this%invp,          &
                      this%rhs,this%nrhs,this%iparm,this%dparm)

  sol = this%rhs

end subroutine solve_pastix_with_rhs

subroutine delete_pastix(this)

  type(pastix_solver) :: this
  sll_int32           :: comm    

  comm = sll_v_world_collective%comm

  ! Call PaStiX clean
  this%iparm(IPARM_START_TASK)       = API_TASK_CLEAN
  this%iparm(IPARM_END_TASK)         = API_TASK_CLEAN
   
  call pastix_fortran(this%pastix_data ,comm,                  &
                      this%ncols,this%colptr,this%row,         &
                      this%avals,this%perm,this%invp,          &
                      this%rhs,this%nrhs,this%iparm,this%dparm)
   
  deallocate(this%colptr)
  deallocate(this%row)
  deallocate(this%avals)
  deallocate(this%perm)
  deallocate(this%invp)
  deallocate(this%rhs)
   
end subroutine delete_pastix

end module sll_m_pastix
