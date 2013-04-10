subroutine mgdrtrsf(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,rc,  &
                    sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,rf,  &
                    comm3dp,myid,neighbor,bd,planetype,IOUT)

use mpi
implicit none 
integer :: sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc
integer :: sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf
integer :: comm3dp,myid,neighbor(26),bd(26),planetype(3),IOUT
real(8) :: rc(sxc-1:exc+1,syc-1:eyc+1,szc-1:ezc+1)
real(8) :: rf(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
!------------------------------------------------------------------------
! Transfer values of the density from a finer to a coarser grid level.
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : gxch1pla, MPI_WAITALL (non-blocking version)
!------------------------------------------------------------------------
integer :: i,j,k,ic,jc,kc,iinc,jinc,kinc,i1,i2,j1,j2,k1,k2
integer :: ireq,req(52)

integer :: status(MPI_STATUS_SIZE,52),ierr

if (nxc.lt.nxf) then
  iinc=2
  i1=1
  i2=0
else 
  iinc=1 
  i1=0
  i2=1
end if
if (nyc.lt.nyf) then
  jinc=2
  j1=1
  j2=0
else
  jinc=1
  j1=0
  j2=1
end if
if (nzc.lt.nzf) then
  kinc=2
  k1=1
  k2=0
else
  kinc=1
  k1=0
  k2=1
end if
do kc=szc,ezc
  k=k1*(2*kc-1)+k2*kc
  do jc=syc,eyc
    j=j1*(2*jc-1)+j2*jc
    do ic=sxc,exc
      i=i1*(2*ic-1)+i2*ic
      rc(ic,jc,kc)=rf(i,j,k)
    end do
  end do
end do
!
! exchange the plane boundary values at coarse level
!

ireq=0

call gxch1pla(sxc,exc,syc,eyc,szc,ezc,rc,comm3dp,neighbor, &
              bd,planetype,req,ireq,IOUT)

call MPI_WAITALL(ireq,req,status,ierr)

return
end
