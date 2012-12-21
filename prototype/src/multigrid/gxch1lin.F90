subroutine gxch1lin(sx,ex,sy,ey,sz,ez,a,comm3d,neighbor, &
                    bd,linetype,req,ireq,IOUT)

use mpi
implicit none 
integer :: sx,ex,sy,ey,sz,ez,comm3d
integer :: neighbor(26),bd(26),linetype(3),req(52),ireq,IOUT
real(8) :: a(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
!------------------------------------------------------------------------
! Subroutine to exchange lines of thickness 1 of boundary data between
! 'myid' and the 12 processes in the directions (i+,j+), (i+,j-),
! (i-,j+), (i-,j-), (i+,k+), (i+,k-), (i-,k+), (i-,k-), (j+,k+),
! (j+,k-), (j-,k+), (j-,k-). Can be used to communicate scalar 
! (p, r,...) as well as vector quantities ((u,v,w), (ut,vt,wt)). 
! The appropriate datatypes have been defined in the 'type_mpi' 
! subroutine.
!
! plane i
!      6|           |8
!      ---------------                k
!       |           |               ^ 
!       |           |               |
!       |    myid   |               |
!       |           |              ------>
!       |           |               |     j
!      ---------------
!      4|           |2
!
! plane i-1                  plane i+1
!       |    16     |              |     25    | 
!      ---------------            --------------- 
!       |           |              |           |
!       |           |              |           |
!     14|           |10          23|           |19
!       |           |              |           |
!       |           |              |           |
!      ---------------            ---------------
!       |    12     |              |     21    |
!
! Code      : tmgd3, test program for 3D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : setdens, gbackin, gdrbsor, gprbsor, gveloc, mgdrestr, 
!             mgdsolver
! Calls     : MPI_ISEND, MPI_IRECV (non-blocking version)
!             MPI_SENDRECV (blocking version)
!------------------------------------------------------------------------
integer :: ierr
!--------------------------non-blocking----------------------------------
! lines (i=constant,j=constant)
! 
! send to 19
!
if (bd(19).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,ey,sz),1,linetype(1),neighbor(19), &
                 6,comm3d,req(ireq),ierr)
end if
!
! receive from 14
!
if (bd(14).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,sy-1,sz),1,linetype(1),neighbor(14), &
                 6,comm3d,req(ireq),ierr)
end if
!
! send to 23
!
if (bd(23).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,sy,sz),1,linetype(1),neighbor(23), &
                 7,comm3d,req(ireq),ierr)
end if
!
! receive from 10
!
if (bd(10).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,ey+1,sz),1,linetype(1),neighbor(10), &
                 7,comm3d,req(ireq),ierr)
end if
!
! send to 14
!
if (bd(14).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,sz),1,linetype(1),neighbor(14), &
                 8,comm3d,req(ireq),ierr)
end if
!
! receive from 19
!
if (bd(19).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,ey+1,sz),1,linetype(1),neighbor(19), &
                 8,comm3d,req(ireq),ierr)
end if
!
! send to 10
!
if (bd(10).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,ey,sz),1,linetype(1),neighbor(10), &
                 9,comm3d,req(ireq),ierr)
end if
!
! receive from 23
!
if (bd(23).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,sy-1,sz),1,linetype(1),neighbor(23), &
                 9,comm3d,req(ireq),ierr)
end if
!------------------------------------------------------------------------
! lines (i=constant,k=constant)
! 
! send to 25
!
if (bd(25).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,sy,ez),1,linetype(2),neighbor(25), &
                 10,comm3d,req(ireq),ierr)
end if
!
! receive from 12
!
if (bd(12).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,sy,sz-1),1,linetype(2),neighbor(12), &
                 10,comm3d,req(ireq),ierr)
end if
!
! send to 21
!
if (bd(21).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,sy,sz),1,linetype(2),neighbor(21), &
                 11,comm3d,req(ireq),ierr)
end if
!
! receive from 16
!
if (bd(16).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,sy,ez+1),1,linetype(2),neighbor(16), &
                 11,comm3d,req(ireq),ierr)
end if
!
! send to 12
!
if (bd(12).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,sz),1,linetype(2),neighbor(12), &
                 12,comm3d,req(ireq),ierr)
end if
!
! receive from 25
!
if (bd(25).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,sy,ez+1),1,linetype(2),neighbor(25), &
                 12,comm3d,req(ireq),ierr)
end if
!
! send to 16
!
if (bd(16).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,ez),1,linetype(2),neighbor(16), &
                 13,comm3d,req(ireq),ierr)
end if
!
! receive from 21
!
if (bd(21).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,sy,sz-1),1,linetype(2),neighbor(21), &
                 13,comm3d,req(ireq),ierr)
end if
!------------------------------------------------------------------------
! lines (j=constant,k=constant)
! 
! send to 8
!
if (bd(8).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,ey,ez),1,linetype(3),neighbor(8), &
                 14,comm3d,req(ireq),ierr)
end if
!
! receive from 4
!
if (bd(4).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,sy-1,sz-1),1,linetype(3),neighbor(4), &
                 14,comm3d,req(ireq),ierr)
end if
!
! send to 2
!
if (bd(2).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,ey,sz),1,linetype(3),neighbor(2), &
                 15,comm3d,req(ireq),ierr)
end if
!
! receive from 6
!
if (bd(6).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,sy-1,ez+1),1,linetype(3),neighbor(6), &
                 15,comm3d,req(ireq),ierr)
end if
!
! send to 4
!
if (bd(4).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,sz),1,linetype(3),neighbor(4), &
                 16,comm3d,req(ireq),ierr)
end if
!
! receive from 8
!
if (bd(8).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,ey+1,ez+1),1,linetype(3),neighbor(8), &
                 16,comm3d,req(ireq),ierr)
end if
!
! send to 6
!
if (bd(6).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,ez),1,linetype(3),neighbor(6), &
                 17,comm3d,req(ireq),ierr)
end if
!
! receive from 2
!
if (bd(2).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,ey+1,sz-1),1,linetype(3),neighbor(2), &
                 17,comm3d,req(ireq),ierr)
end if

return
end

subroutine gxch1lin(a,comm2d,sx,ex,sy,ey,neighbor,bd,
                    idatatype,jdatatype,IOUT)
#include "mgd2.h"
include "mpif.h"
integer sx,ex,sy,ey,IOUT
REALN a(sx-1:ex+1,sy-1:ey+1)
integer comm2d,neighbor(8),bd(8),idatatype,jdatatype
!------------------------------------------------------------------------
! Subroutine to exchange one-lines (one row and one column) of 
! boundary data between "directly" neighboring processes. This subroutine 
! can be used to exchange scalar as well as vector variables; it suffices 
! to pass the right argument for the datatypes: 
! idatatype, jdatatype -> r, p, tmp...
! i11datatype, j11datatype -> (u,v)
! i12datatype, j12datatype -> (ut,vt)
!
! If other quantities need to be passed, appropriate datatypes for
! them have to be defined first in type_mpi.
!
! 'neighbor' and 'bd' arrays:
!
!       |     7     | 
!       |           |
!    ------------------
!       |           |
!       |           |
!     5 |   myid    | 1 
!       |           |
!       |           |
!    ------------------
!       |           |
!       |     3     |  
!
! Code      : tmgd2
! Called in : mgdrelax, mgdrelax, mgdrestr, mgdrtrsf
! Calls     : MPI_ISEND, MPI_IRECV, MPI_WAITALL (non-blocking version)
!             MPI_SENDRECV (blocking version)
!------------------------------------------------------------------------
# if NBLOCKGR
integer req(8),status(MPI_STATUS_SIZE,8),ireq,ierr
# else
integer status(MPI_STATUS_SIZE),ierr
# endif
# if cdebug
# if NBLOCKGR
integer nc
# endif
double precision tinitial
tinitial=MPI_WTIME()
# endif
# if NBLOCKGR
!--------------------------non-blocking----------------------------------
ireq=0
!
! send to 1
!
if (bd(1).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,sy),1,jdatatype,neighbor(1), &
                 0,comm2d,req(ireq),ierr)
end if
!
! receive from 5
!
if (bd(5).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,sy),1,jdatatype,neighbor(5), &
                 0,comm2d,req(ireq),ierr)
end if
!
! send to 3
!
if (bd(3).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy),1,idatatype,neighbor(3), &
                 1,comm2d,req(ireq),ierr)
end if
!
! receive from 7
! 
if (bd(7).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,ey+1),1,idatatype,neighbor(7), &
                 1,comm2d,req(ireq),ierr)
end if
!
! send to 5
!
if (bd(5).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy),1,jdatatype,neighbor(5), &
                 1,comm2d,req(ireq),ierr)
end if
!
! receive from 1
!
if (bd(1).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,sy),1,jdatatype,neighbor(1), &
                 1,comm2d,req(ireq),ierr)
end if
!
! send to 7
!
if (bd(7).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,ey),1,idatatype,neighbor(7), &
                 0,comm2d,req(ireq),ierr)
end if
!
! receive from 3
!
if (bd(3).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,sy-1),1,idatatype,neighbor(3), &
                 0,comm2d,req(ireq),ierr)
end if
!
! wait for all the messages to be sent and received before going on.
!
call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
nc=4-(bd(1)+bd(3)+bd(5)+bd(7))
nisend(1,1)=nisend(1,1)+nc
nirecv(1,1)=nirecv(1,1)+nc
nwaitall=nwaitall+1
# endif
# else
!----------------------------blocking------------------------------------
! send to 1 and receive from 5
!
call MPI_SENDRECV(a(ex,sy),1,jdatatype,neighbor(1),0, &
                  a(sx-1,sy),1,jdatatype,neighbor(5),0, &
                  comm2d,status,ierr)
!
! send to 3 and receive from 7
!
call MPI_SENDRECV(a(sx,sy),1,idatatype,neighbor(3),1, &
                  a(sx,ey+1),1,idatatype,neighbor(7),1, &
                  comm2d,status,ierr)
!
! send to 5 and receive from 1
!
call MPI_SENDRECV(a(sx,sy),1,jdatatype,neighbor(5),1, &
                  a(ex+1,sy),1,jdatatype,neighbor(1),1, &
                  comm2d,status,ierr)
!
! send to 7 and receive from 3
!
call MPI_SENDRECV(a(sx,ey),1,idatatype,neighbor(7),0, &
                  a(sx,sy-1),1,idatatype,neighbor(3),0, &
                  comm2d,status,ierr)
# if cdebug
nsendrecv(1,1)=nsendrecv(1,1)+4
# endif
# endif
c
# if cdebug
timing(59)=timing(59)+MPI_WTIME()-tinitial
# endif
return
end
