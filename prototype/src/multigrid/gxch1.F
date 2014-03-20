      subroutine gxch1cor(a,comm2d,sx,ex,sy,ey,neighbor,bd,
     1                    ijdatatype,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1)
      integer comm2d,neighbor(8),bd(8),ijdatatype
      integer ierr
c------------------------------------------------------------------------
c Subroutine to exchange one corner point of boundary data between 
c "diagonally" neighboring processes. This subroutineccan be used to 
c exchange scalar as well as vector variables; it suffices to pass the 
c right argument for the datatypes:
c ijdatatype -> r, p, tmp...
c ij11datatype -> (u,v)
c ij12datatype -> (ut,vt)
c
c 'neighbor' and 'bd' arrays:
c
c     6 |           | 8
c       |           |
c    ------------------
c       |           |
c       |           |
c       |   myid    |   
c       |           |
c       |           |
c    ------------------
c       |           |
c     4 |           | 2
c
c Code      : tmgd2
c Called in : mgdrestr, mgdsolver
c Calls     : MPI_ISEND, MPI_IRECV, MPI_WAITALL (non-blocking version)
c             MPI_SENDRECV (blocking version)
c----------------------------------------------------------------------- 
# if NBLOCKGR
      integer req(8),status(MPI_STATUS_SIZE,8),ireq
# else
      integer status(MPI_STATUS_SIZE)
# endif
# if cdebug
# if NBLOCKGR
      integer nc
# endif
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
# if NBLOCKGR
c--------------------------non-blocking----------------------------------
      ireq=0
c
c send to 2
c
      if (bd(2).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(ex,sy),1,ijdatatype,neighbor(2),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c receive from 6
c
      if (bd(6).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx-1,ey+1),1,ijdatatype,neighbor(6),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c send to 4
c
      if (bd(4).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,sy),1,ijdatatype,neighbor(4),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c receive from 8
c
      if (bd(8).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(ex+1,ey+1),1,ijdatatype,neighbor(8),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c send to 6
c
      if (bd(6).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,ey),1,ijdatatype,neighbor(6),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c receive from 2
c
      if (bd(2).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(ex+1,sy-1),1,ijdatatype,neighbor(2),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c send to 8
c
      if (bd(8).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(ex,ey),1,ijdatatype,neighbor(8),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c receive from 4
c
      if (bd(4).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx-1,sy-1),1,ijdatatype,neighbor(4),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c wait for all the messages to be sent and received before going on.
c
      call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
      nc=4-(bd(2)+bd(4)+bd(6)+bd(8))
      nisend(2,1)=nisend(2,1)+nc
      nirecv(2,1)=nirecv(2,1)+nc
      nwaitall=nwaitall+1
# endif
# else
c----------------------------blocking------------------------------------
c send to 2 and receive from 6
c
      call MPI_SENDRECV(a(ex,sy),1,ijdatatype,neighbor(2),0,
     1                  a(sx-1,ey+1),1,ijdatatype,neighbor(6),0,
     2                  comm2d,status,ierr)
c
c send to 4 and receive from 8
c
      call MPI_SENDRECV(a(sx,sy),1,ijdatatype,neighbor(4),1,
     1                  a(ex+1,ey+1),1,ijdatatype,neighbor(8),1,
     2                  comm2d,status,ierr)
c
c send to 6 and receive from 2
c
      call MPI_SENDRECV(a(sx,ey),1,ijdatatype,neighbor(6),1,
     1                  a(ex+1,sy-1),1,ijdatatype,neighbor(2),1,
     2                  comm2d,status,ierr)
c
c send to 8 and receive from 4
c
      call MPI_SENDRECV(a(ex,ey),1,ijdatatype,neighbor(8),0,
     1                  a(sx-1,sy-1),1,ijdatatype,neighbor(4),0,
     2                  comm2d,status,ierr)
# if cdebug
      nsendrecv(2,1)=nsendrecv(2,1)+4
# endif
# endif
c
# if cdebug
      timing(60)=timing(60)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine gxch1lin(a,comm2d,sx,ex,sy,ey,neighbor,bd,
     1                    idatatype,jdatatype,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1)
      integer comm2d,neighbor(8),bd(8),idatatype,jdatatype
      integer ierr
c------------------------------------------------------------------------
c Subroutine to exchange one-lines (one row and one column) of 
c boundary data between "directly" neighboring processes. This subroutine 
c can be used to exchange scalar as well as vector variables; it suffices 
c to pass the right argument for the datatypes: 
c idatatype, jdatatype -> r, p, tmp...
c i11datatype, j11datatype -> (u,v)
c i12datatype, j12datatype -> (ut,vt)
c
c If other quantities need to be passed, appropriate datatypes for
c them have to be defined first in type_mpi.
c
c 'neighbor' and 'bd' arrays:
c
c       |     7     | 
c       |           |
c    ------------------
c       |           |
c       |           |
c     5 |   myid    | 1 
c       |           |
c       |           |
c    ------------------
c       |           |
c       |     3     |  
c
c Code      : tmgd2
c Called in : mgdrelax, mgdrelax, mgdrestr, mgdrtrsf
c Calls     : MPI_ISEND, MPI_IRECV, MPI_WAITALL (non-blocking version)
c             MPI_SENDRECV (blocking version)
c------------------------------------------------------------------------
# if NBLOCKGR
      integer req(8),status(MPI_STATUS_SIZE,8),ireq
# else
      integer status(MPI_STATUS_SIZE)
# endif
# if cdebug
# if NBLOCKGR
      integer nc
# endif
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
# if NBLOCKGR
c--------------------------non-blocking----------------------------------
      ireq=0
c
c send to 1
c
      if (bd(1).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(ex,sy),1,jdatatype,neighbor(1),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c receive from 5
c
      if (bd(5).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx-1,sy),1,jdatatype,neighbor(5),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c send to 3
c
      if (bd(3).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,sy),1,idatatype,neighbor(3),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c receive from 7
c 
      if (bd(7).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx,ey+1),1,idatatype,neighbor(7),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c send to 5
c
      if (bd(5).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,sy),1,jdatatype,neighbor(5),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c receive from 1
c
      if (bd(1).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(ex+1,sy),1,jdatatype,neighbor(1),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c send to 7
c
      if (bd(7).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,ey),1,idatatype,neighbor(7),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c receive from 3
c
      if (bd(3).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx,sy-1),1,idatatype,neighbor(3),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c wait for all the messages to be sent and received before going on.
c
      call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
      nc=4-(bd(1)+bd(3)+bd(5)+bd(7))
      nisend(1,1)=nisend(1,1)+nc
      nirecv(1,1)=nirecv(1,1)+nc
      nwaitall=nwaitall+1
# endif
# else
c----------------------------blocking------------------------------------
c send to 1 and receive from 5
c
      call MPI_SENDRECV(a(ex,sy),1,jdatatype,neighbor(1),0,
     1                  a(sx-1,sy),1,jdatatype,neighbor(5),0,
     2                  comm2d,status,ierr)
c
c send to 3 and receive from 7
c
      call MPI_SENDRECV(a(sx,sy),1,idatatype,neighbor(3),1,
     1                  a(sx,ey+1),1,idatatype,neighbor(7),1,
     2                  comm2d,status,ierr)
c
c send to 5 and receive from 1
c
      call MPI_SENDRECV(a(sx,sy),1,jdatatype,neighbor(5),1,
     1                  a(ex+1,sy),1,jdatatype,neighbor(1),1,
     2                  comm2d,status,ierr)
c
c send to 7 and receive from 3
c
      call MPI_SENDRECV(a(sx,ey),1,idatatype,neighbor(7),0,
     1                  a(sx,sy-1),1,idatatype,neighbor(3),0,
     2                  comm2d,status,ierr)
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

