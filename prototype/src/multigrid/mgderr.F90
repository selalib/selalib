module mgderr
#include "sll_working_precision.h"
use mpi
implicit none 

contains

subroutine mgderr_3d(relmax,sxm,exm,sym,eym,szm,ezm,phio,phin,comm3d)
integer :: sxm,exm,sym,eym,szm,ezm,comm3d
real(8) :: relmax
real(8) :: phio(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1)
real(8) :: phin(sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1)
!------------------------------------------------------------------------
! Calculate the error between the new and old iterates of phi and 
! save the new iterate into the phio array.
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : MPI_ALLREDUCE
!------------------------------------------------------------------------
real(8) :: locval(2),totval(2)
integer :: i,j,k,ierr
!
! calculate local error
!
locval(1)=0.0d0
locval(2)=0.0d0
do k=szm,ezm
  do j=sym,eym
    do i=sxm,exm
      locval(1)=max(locval(1),abs(phin(i,j,k)))
      locval(2)=max(locval(2),abs(phin(i,j,k)-phio(i,j,k)))
    end do
  end do
end do
!
! global reduce across all processes
!
call MPI_ALLREDUCE(locval,totval,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)

if (totval(1).gt.0.0) then
  relmax=totval(2)/totval(1)
else
  relmax=0.0d0
end if
!
! save new values into ouput array
!
do k=szm-1,ezm+1
  do j=sym-1,eym+1
    do i=sxm-1,exm+1
      phio(i,j,k)=phin(i,j,k)
    end do
  end do
end do

return
end subroutine

subroutine mgderr_2d(relmax,sxm,exm,sym,eym,phio,phin,comm2d)
# include "mgd2.h"
integer sxm,exm,sym,eym,comm2d
REALN relmax
REALN phio(sxm-1:exm+1,sym-1:eym+1),phin(sxm-1:exm+1,sym-1:eym+1)
!------------------------------------------------------------------------
! Calculate the error between the new and old iterates of phi and 
! save the new iterate into the phio array.
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : MPI_ALLREDUCE
!------------------------------------------------------------------------
REALN phloc,reloc
integer i,j,ierr
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
# if double_precision
call MPI_ALLREDUCE(reloc,relmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm2d,ierr)
# else
call MPI_ALLREDUCE(reloc,relmax,1,MPI_REAL,MPI_MAX,comm2d,ierr)
# endif
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

end module mgderr
