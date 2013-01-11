subroutine gxch1cor(a,comm2d,sx,ex,sy,ey,neighbor,bd,ijdatatype)
use mpi
#include "sll_working_precision.h"
#include "mgd2.h"
sll_int32  :: sx,ex,sy,ey
sll_real64 :: a(sx-1:ex+1,sy-1:ey+1)
sll_int32  :: comm2d,neighbor(8),bd(8),ijdatatype
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
integer :: req(8),status(MPI_STATUS_SIZE,8),ireq,ierr
# else
integer :: status(MPI_STATUS_SIZE),ierr
# endif
# if DEBUG
# if NBLOCKGR
integer :: nc
# endif
double precision :: tinitial
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
# if DEBUG
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
# if DEBUG
nsendrecv(2,1)=nsendrecv(2,1)+4
# endif
# endif

# if DEBUG
timing(60)=timing(60)+MPI_WTIME()-tinitial
# endif
return
end subroutine
