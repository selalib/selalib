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


subroutine gxch1cor(a,comm2d,sx,ex,sy,ey,neighbor,bd,ijdatatype,IOUT)
# include "mgd2.h"
include "mpif.h"
integer sx,ex,sy,ey,IOUT
REALN a(sx-1:ex+1,sy-1:ey+1)
integer comm2d,neighbor(8),bd(8),ijdatatype
!------------------------------------------------------------------------
! Subroutine to exchange one corner point of boundary data between 
! "diagonally" neighboring processes. This subroutineccan be used to 
! exchange scalar as well as vector variables; it suffices to pass the 
! right argument for the datatypes:
! ijdatatype -> r, p, tmp...
! ij11datatype -> (u,v)
! ij12datatype -> (ut,vt)
!
! 'neighbor' and 'bd' arrays:
!
!     6 |           | 8
!       |           |
!    ------------------
!       |           |
!       |           |
!       |   myid    |   
!       |           |
!       |           |
!    ------------------
!       |           |
!     4 |           | 2
!
! Code      : tmgd2
! Called in : mgdrestr, mgdsolver
! Calls     : MPI_ISEND, MPI_IRECV, MPI_WAITALL (non-blocking version)
!             MPI_SENDRECV (blocking version)
!----------------------------------------------------------------------- 
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
! send to 2
!
if (bd(2).eq.0) then
   ireq=ireq+1
   call MPI_ISEND(a(ex,sy),1,ijdatatype,neighbor(2), &
                  0,comm2d,req(ireq),ierr)
end if
!
! receive from 6
!
if (bd(6).eq.0) then
   ireq=ireq+1
   call MPI_IRECV(a(sx-1,ey+1),1,ijdatatype,neighbor(6), &
                  0,comm2d,req(ireq),ierr)
end if
!
! send to 4
!
if (bd(4).eq.0) then
   ireq=ireq+1
   call MPI_ISEND(a(sx,sy),1,ijdatatype,neighbor(4), &
                  1,comm2d,req(ireq),ierr)
end if
!
! receive from 8
!
if (bd(8).eq.0) then
   ireq=ireq+1
   call MPI_IRECV(a(ex+1,ey+1),1,ijdatatype,neighbor(8), &
                    1,comm2d,req(ireq),ierr)
end if
!
! send to 6
!
if (bd(6).eq.0) then
   ireq=ireq+1
   call MPI_ISEND(a(sx,ey),1,ijdatatype,neighbor(6), &
                  1,comm2d,req(ireq),ierr)
end if
!
! receive from 2
!
if (bd(2).eq.0) then
   ireq=ireq+1
   call MPI_IRECV(a(ex+1,sy-1),1,ijdatatype,neighbor(2), &
                  1,comm2d,req(ireq),ierr)
end if
!
! send to 8
!
if (bd(8).eq.0) then
   ireq=ireq+1
   call MPI_ISEND(a(ex,ey),1,ijdatatype,neighbor(8), &
                  0,comm2d,req(ireq),ierr)
end if
!
! receive from 4
!
if (bd(4).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,sy-1),1,ijdatatype,neighbor(4), &
                      0,comm2d,req(ireq),ierr)
end if
!
! wait for all the messages to be sent and received before going on.
!
call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
nc=4-(bd(2)+bd(4)+bd(6)+bd(8))
nisend(2,1)=nisend(2,1)+nc
nirecv(2,1)=nirecv(2,1)+nc
nwaitall=nwaitall+1
# endif
# else
!----------------------------blocking------------------------------------
! send to 2 and receive from 6
!
call MPI_SENDRECV(a(ex,sy),1,ijdatatype,neighbor(2),0,      &
                  a(sx-1,ey+1),1,ijdatatype,neighbor(6),0,  &
                  comm2d,status,ierr)
!
! send to 4 and receive from 8
!
call MPI_SENDRECV(a(sx,sy),1,ijdatatype,neighbor(4),1,      &
                  a(ex+1,ey+1),1,ijdatatype,neighbor(8),1,  &
                  comm2d,status,ierr)
!
! send to 6 and receive from 2
!
call MPI_SENDRECV(a(sx,ey),1,ijdatatype,neighbor(6),1,      &
                  a(ex+1,sy-1),1,ijdatatype,neighbor(2),1,  &
                  comm2d,status,ierr)
!
! send to 8 and receive from 4
!
call MPI_SENDRECV(a(ex,ey),1,ijdatatype,neighbor(8),0,      &
                  a(sx-1,sy-1),1,ijdatatype,neighbor(4),0,  &
                  comm2d,status,ierr)
# if cdebug
nsendrecv(2,1)=nsendrecv(2,1)+4
# endif
# endif

# if cdebug
timing(60)=timing(60)+MPI_WTIME()-tinitial
# endif
return
end
