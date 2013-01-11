subroutine gxch1cor(sx,ex,sy,ey,sz,ez,a,comm3d,neighbor, &
                    bd,cornertype,req,ireq,IOUT)
use mpi
implicit none
#include "mgd3.h"
integer :: sx,ex,sy,ey,sz,ez,comm3d
integer :: neighbor(26),bd(26),cornertype,req(52),ireq,IOUT
real(8) :: a(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
!------------------------------------------------------------------------
! Subroutine to exchange 1*1*1 corners (i.e. points) of boundary data 
! between 'myid' and the 8 processes (i+,j+,k+), (i+,j+,k-), 
! (i+,j-,k+), (i+,j-,k-), (i-,j+,k+), (i-,j+,k-), (i-,j-,k+),
! (i-,j-,k-). Can be used to communicate scalar (p, r,...) as well as 
! vector quantities ((u,v,w), (ut,vt,wt)). The appropriate datatypes 
! have been defined in the 'type_mpi' subroutine.
!
!
! plane i
!       |           | 
!      ---------------                k
!       |           |               ^ 
!       |           |               |
!       |    myid   |               |
!       |           |              ------>
!       |           |               |     j
!      ---------------
!       |           | 
!
! plane i-1                  plane i+1
!     15|           |17          24|           |26 
!      ---------------            --------------- 
!       |           |              |           |
!       |           |              |           |
!       |           |              |           |
!       |           |              |           |
!       |           |              |           |
!      ---------------            ---------------
!     13|           |11          22|           |20
!
! Code      : tmgd3, test program for 3D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : gbackin, setdens, gdrbsor, gprbsor, gveloc, mgdrestr, 
!             mgdsolver
! Calls     : MPI_ISEND, MPI_IRECV (non-blocking version)
!             MPI_SENDRECV (blocking version)
!------------------------------------------------------------------------
integer ierr
!--------------------------non-blocking----------------------------------
!
! send to 26
!
if (bd(26).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,ey,ez),1,cornertype,neighbor(26), &
  &              18,comm3d,req(ireq),ierr)
end if
!
! receive from 13
!
if (bd(13).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,sy-1,sz-1),1,cornertype,neighbor(13), &
  &              18,comm3d,req(ireq),ierr)
end if
!
! send to 20
!
if (bd(20).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,ey,sz),1,cornertype,neighbor(20), &
  &              19,comm3d,req(ireq),ierr)
end if
!
! receive from 15
!
if (bd(15).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,sy-1,ez+1),1,cornertype,neighbor(15), &
  &              19,comm3d,req(ireq),ierr)
end if
!
! send to 24
!
if (bd(24).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,sy,ez),1,cornertype,neighbor(24), &
  &              20,comm3d,req(ireq),ierr)
end if
!
! receive from 11
!
if (bd(11).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,ey+1,sz-1),1,cornertype,neighbor(11), &
  &              20,comm3d,req(ireq),ierr)
end if
!
! send to 22
!
if (bd(22).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,sy,sz),1,cornertype,neighbor(22), &
  &              21,comm3d,req(ireq),ierr)
end if
!
! receive from 17
!
if (bd(17).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,ey+1,ez+1),1,cornertype,neighbor(17), &
  &              21,comm3d,req(ireq),ierr)
end if
!
! send to 13
!
if (bd(13).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,sz),1,cornertype,neighbor(13), &
  &              22,comm3d,req(ireq),ierr)
end if
!
! receive from 26
!
if (bd(26).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,ey+1,ez+1),1,cornertype,neighbor(26), &
  &              22,comm3d,req(ireq),ierr)
end if
!
! send to 15
!
if (bd(15).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,ez),1,cornertype,neighbor(15), &
  &              23,comm3d,req(ireq),ierr)
end if
!
! receive from 20
!
if (bd(20).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,ey+1,sz-1),1,cornertype,neighbor(20), &
  &              23,comm3d,req(ireq),ierr)
end if
!
! send to 11
!
if (bd(11).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,ey,sz),1,cornertype,neighbor(11), &
  &              24,comm3d,req(ireq),ierr)
end if
!
! receive from 24
!
if (bd(24).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,sy-1,ez+1),1,cornertype,neighbor(24), &
  &              24,comm3d,req(ireq),ierr)
end if
!
! send to 17
!
if (bd(17).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,ey,ez),1,cornertype,neighbor(17), &
  &              25,comm3d,req(ireq),ierr)
end if
!
! receive from 22
!
if (bd(22).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,sy-1,sz-1),1,cornertype,neighbor(22), &
  &              25,comm3d,req(ireq),ierr)
end if
return
end
