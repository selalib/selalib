module gxch1_2d
use mpi
#include "sll_working_precision.h"

# if NBLOCKGR
sll_int32 :: nisend(2,3)
sll_int32 :: nirecv(2,3)
sll_int32 :: nwait
sll_int32 :: nwaitall
# else
sll_int32  :: nsendrecv(2,3)
# endif

sll_int32  :: nallreduce
sll_int32  :: nalltoall
sll_int32  :: nreduce

sll_int32  :: nisendfr
sll_int32  :: nirecvfr
sll_int32  :: nwaitallfr
sll_int32  :: nsteptiming
logical    :: nocterr

contains

subroutine gxch1cor(a,comm2d,sx,ex,sy,ey,neighbor,bd,ijdatatype)
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
!timing(60)=timing(60)+MPI_WTIME()-tinitial
# endif
return
end subroutine


subroutine gxch1lin(a,comm2d,sx,ex,sy,ey,neighbor,bd, &
                    idatatype,jdatatype)
sll_int32  :: sx,ex,sy,ey
sll_real64 :: a(sx-1:ex+1,sy-1:ey+1)
sll_int32  :: comm2d,neighbor(8),bd(8),idatatype,jdatatype
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
sll_int32 :: req(8),status(MPI_STATUS_SIZE,8),ireq,ierr
# else
sll_int32 :: status(MPI_STATUS_SIZE),ierr
# endif
# if DEBUG
# if NBLOCKGR
sll_int32 :: nc
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
# if DEBUG
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
# if DEBUG
nsendrecv(1,1)=nsendrecv(1,1)+4
# endif
# endif

# if DEBUG
!timing(59)=timing(59)+MPI_WTIME()-tinitial
# endif

end subroutine

end module gxch1_2d

