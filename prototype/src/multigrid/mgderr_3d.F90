subroutine mgderr(relmax,sxm,exm,sym,eym,szm,ezm,phio,phin,comm3d,IOUT)

use mpi
implicit none 

integer :: sxm,exm,sym,eym,szm,ezm,comm3d,IOUT
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
call MPI_ALLREDUCE(locval,totval,2,MPI_REAL8,MPI_MAX,comm3d,ierr)

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
end
