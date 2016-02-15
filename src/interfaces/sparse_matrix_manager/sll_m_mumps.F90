module sll_m_mumps
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define MPI_MASTER 0
#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_m_collective, only:  &
  sll_s_boot_collective,     &
  sll_f_get_collective_rank, &
  sll_f_get_collective_size, &
  sll_v_world_collective

implicit none


public :: &
  mumps_solver, &
  initialize, &
  factorize, &
  solve, &
  delete

private

INCLUDE 'dmumps_struc.h'

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
type :: mumps_solver
 type(dmumps_struc) :: mumps_par !< MUMPS structure 
end type mumps_solver

interface initialize
   module procedure init_mumps
end interface initialize

interface solve
   module procedure solve_mumps_with_rhs
end interface solve

interface factorize
   module procedure factorize_mumps
end interface factorize

interface delete
   module procedure free_mumps
end interface delete

contains

subroutine init_mumps(self,n,nnzeros,row_ptr,col_ind,val)

  type(mumps_solver)               :: self
  sll_int32,    intent(in)         :: n
  sll_int32,    intent(in)         :: nnzeros
  sll_int32,    intent(in), target :: row_ptr(:)
  sll_int32,    intent(in), target :: col_ind(:)
  sll_real64,   intent(in), target :: val(:)
  sll_int32                        :: error
  sll_int32                        :: prank    
  sll_int32                        :: psize

  if( .not. associated(sll_v_world_collective)) then
     call sll_s_boot_collective()
  end if
  prank = sll_f_get_collective_rank( sll_v_world_collective )
  psize = sll_f_get_collective_size( sll_v_world_collective )
  ! Define a communicator for the package.
  self%mumps_par%comm = sll_v_world_collective%comm

  !  Initialize an instance of the package
  !  for L U factorization (sym = 0, with working host)
  self%mumps_par%JOB = -1
  self%mumps_par%SYM = 0
  self%mumps_par%PAR = 1
  call dmumps(self%mumps_par)

  if (self%mumps_par%INFOG(1).LT.0) then
    write(6,'(A,A,I6,A,I9)') " ERROR RETURN: ", &
              "  mumps_par%INFOG(1)= ", self%mumps_par%INFOG(1),  &
              "  mumps_par%INFOG(2)= ", self%mumps_par%INFOG(2) 
    stop
  end if

  self%mumps_par%N  = n
  self%mumps_par%NZ = nnzeros
  allocate( self%mumps_par%IRN ( self%mumps_par%NZ ) )
  allocate( self%mumps_par%JCN ( self%mumps_par%NZ ) )
  allocate( self%mumps_par%A(    self%mumps_par%NZ ) )
  allocate( self%mumps_par%RHS ( self%mumps_par%N  ) )
  !mumps_par%IRN
  !mumps_par%JCN
  !mumps_par%A(I)
  !mumps_par%RHS(I)

end subroutine init_mumps

subroutine factorize_mumps(self)

  type(mumps_solver) :: self

end subroutine factorize_mumps

subroutine solve_mumps_with_rhs(self, rhs, sol)

  type(mumps_solver)      :: self
  sll_real64, dimension(:) :: rhs
  sll_real64, dimension(:) :: sol

end subroutine solve_mumps_with_rhs

subroutine free_mumps(self)

  type(mumps_solver) :: self

end subroutine free_mumps

end module sll_m_mumps
