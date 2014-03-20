      subroutine gxch1cor(a,comm2d,sx,ex,sy,ey,neighbor,bd,
     1                    ijdatatype,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1)
      integer comm2d,neighbor(8),bd(8),ijdatatype
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

