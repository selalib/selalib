subroutine mgderr(relmax,sxm,exm,sym,eym,phio,phin,comm2d)
use mpi
#include "sll_working_precision.h"
#include "mgd2.h"
sll_int32  :: sxm,exm,sym,eym,comm2d
sll_real64 :: relmax
sll_real64 :: phio(sxm-1:exm+1,sym-1:eym+1),phin(sxm-1:exm+1,sym-1:eym+1)
!------------------------------------------------------------------------
! Calculate the error between the new and old iterates of phi and 
! save the new iterate into the phio array.
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : MPI_ALLREDUCE
!------------------------------------------------------------------------
sll_real64 :: phloc,reloc
sll_int32  :: i,j,ierr
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif
!
! calculate local error
!
phloc=0.0d0
reloc=0.0d0
do j=sym,eym
  do i=sxm,exm
    phloc=max(phloc,abs(phin(i,j)))
    reloc=max(reloc,abs(phin(i,j)-phio(i,j)))
  end do
end do
if (phloc.gt.0.0) then
  reloc=reloc/phloc
else
  reloc=0.0d0
end if
!
! global reduce across all processes
!
call MPI_ALLREDUCE(reloc,relmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d,ierr)
# if cdebug
nallreduce=nallreduce+1
# endif

! save new values into ouput array

do j=sym-1,eym+1
  do i=sxm-1,exm+1
    phio(i,j)=phin(i,j)
  end do
end do

# if cdebug
timing(94)=timing(94)+MPI_WTIME()-tinitial
# endif

end subroutine
