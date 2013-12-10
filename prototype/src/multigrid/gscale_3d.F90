subroutine gscale(sx,ex,sy,ey,sz,ez,a,avo,acorr,comm3d,nx,ny,nz)

use mpi
implicit none 
integer :: sx,ex,sy,ey,sz,ez,nx,ny,nz
real(8) :: a(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),avo,acorr
integer :: comm3d
!------------------------------------------------------------------------
! Rescale the field a so that its average inside the domain
! remains constant and equal to avo. For the density,avo should
! be rro, this ensures conservation of mass. For the pressure,
! avo should be 0 so that the average pressure does not drift
! away from 0, which is the initial value. For the density, also
! eliminate the numerical overshoots and undershoots in the solution,
! so that can consider higher density ratios; this feature is active
! only if the compiler directive UNDERSHOOT is set to 1 in the
! "compdir.inc" file.
!
! Code      : tmgd3, test program for 3D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver, gdrbsor, gprbsor
! Calls     : MPI_ALLREDUCE
!------------------------------------------------------------------------
real(8) :: avloc,av
integer :: i,j,k,ierr

!
! determine average value
!
avloc=0.0d0
do k=sz,ez
  do j=sy,ey
    do i=sx,ex
      avloc=avloc+a(i,j,k)
    end do
  end do
end do
!
! global reduce across all process
!
call MPI_ALLREDUCE(avloc,av,1,MPI_REAL8,MPI_SUM,comm3d,ierr)

av=av/float(nx*ny*nz)
!
! do correction
!
acorr=avo-av
do k=sz,ez
  do j=sy,ey
    do i=sx,ex
      a(i,j,k)=a(i,j,k)+acorr
    end do
  end do
end do

return
end 

