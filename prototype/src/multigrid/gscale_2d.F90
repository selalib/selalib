subroutine gscale(sx,ex,sy,ey,a,avo,acorr,comm2d,nx,ny)
use mpi
#include "sll_working_precision.h"

sll_int32  :: sx,ex,sy,ey,nx,ny
sll_real64 :: a(sx-1:ex+1,sy-1:ey+1),avo,acorr
sll_int32  :: comm2d
!------------------------------------------------------------------------
! Rescale the field a so that its average inside the domain
! remains constant and equal to avo. For the density,avo should
! be rro, this ensures conservation of mass. For the pressure,
! avo should be 0 so that the average pressure does not drift
! away from 0, which is the initial value.
!
! Code      : tmgd2
! Called in : mgdsolver
! Calls     : MPI_ALLREDUCE
!------------------------------------------------------------------------
sll_real64 :: avloc,av
sll_int32  :: i,j,ierr

! determine average value
avloc=0.0d0
do j=sy,ey
   do i=sx,ex
      avloc=avloc+a(i,j)
   end do
end do
!
! global reduce across all process
!
call MPI_ALLREDUCE(avloc,av,1,MPI_REAL8,MPI_SUM,comm2d,ierr)
# if DEBUG
nallreduce=nallreduce+1
# endif
av=av/float(nx*ny)
!
! do correction
!
acorr=avo-av
do j=sy,ey
  do i=sx,ex
    a(i,j)=a(i,j)+acorr
  end do
end do

end subroutine 
