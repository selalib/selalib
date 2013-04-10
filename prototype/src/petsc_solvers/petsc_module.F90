module petsc_module
#include <finclude/petscdmdef.h>
use petscdmdef
!#include "finclude/petscdef.h"
!use petsc
implicit none
private

type, public    :: userctx
   DM           :: da
   PetscInt     :: xs,xe,xm,gxs,gxe,gxm
   PetscInt     :: ys,ye,ym,gys,gye,gym
   PetscInt     :: mx,my
   PetscMPIInt  :: rank
end type userctx

public :: hello_petsc

contains

subroutine hello_petsc()
PetscInt :: prank, psize
PetscErrorCode :: ierr
character(len=40) :: chaine
call MPI_Comm_split(MPI_COMM_WORLD,mod(prank,2),0,PETSC_COMM_WORLD,ierr)
!Every PETSc routine should begin with the PetscInitialize() routine.
call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
call MPI_Comm_size(PETSC_COMM_WORLD,psize,ierr)
call MPI_Comm_rank(PETSC_COMM_WORLD,prank,ierr)
write(chaine,"('Hello World from node',i2)") prank
!     Here we would like to print only one message that represents all
!     the processes in the group. 
write(6,100) psize,prank
call PetscSynchronizedPrintf(PETSC_COMM_WORLD, chaine//'\n', ierr)
!call PetscSynchronizedPrintf(PETSC_COMM_WORLD,chaine//"\n", ierr)
call PetscSynchronizedFlush(PETSC_COMM_WORLD,ierr)
call PetscFinalize(ierr)

 100  format("No of Procs PETSc = ",i4," rank = ",i4)

end subroutine hello_petsc

end module petsc_module
