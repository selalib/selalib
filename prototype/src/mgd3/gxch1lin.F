      subroutine gxch1lin(sx,ex,sy,ey,sz,ez,a,comm3d,neighbor,
     1                    bd,linetype,req,ireq,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,sz,ez,comm3d
      integer neighbor(26),bd(26),linetype(3),req(52),ireq,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
c------------------------------------------------------------------------
c Subroutine to exchange lines of thickness 1 of boundary data between
c 'myid' and the 12 processes in the directions (i+,j+), (i+,j-),
c (i-,j+), (i-,j-), (i+,k+), (i+,k-), (i-,k+), (i-,k-), (j+,k+),
c (j+,k-), (j-,k+), (j-,k-). Can be used to communicate scalar 
c (p, r,...) as well as vector quantities ((u,v,w), (ut,vt,wt)). 
c The appropriate datatypes have been defined in the 'type_mpi' 
c subroutine.
c
c plane i
c      6|           |8
c      ---------------                k
c       |           |               ^ 
c       |           |               |
c       |    myid   |               |
c       |           |              ------>
c       |           |               |     j
c      ---------------
c      4|           |2
c
c plane i-1                  plane i+1
c       |    16     |              |     25    | 
c      ---------------            --------------- 
c       |           |              |           |
c       |           |              |           |
c     14|           |10          23|           |19
c       |           |              |           |
c       |           |              |           |
c      ---------------            ---------------
c       |    12     |              |     21    |
c
c Code      : tmgd3, test program for 3D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : setdens, gbackin, gdrbsor, gprbsor, gveloc, mgdrestr, 
c             mgdsolver
c Calls     : MPI_ISEND, MPI_IRECV (non-blocking version)
c             MPI_SENDRECV (blocking version)
c------------------------------------------------------------------------
      integer ierr
# if !NBLOCKGR
      integer status(MPI_STATUS_SIZE)
# endif
# if cdebug
# if NBLOCKGR
      integer nc
# endif
      double precision tinitial,tmpi
      tinitial=MPI_WTIME()
# endif
# if NBLOCKGR
c--------------------------non-blocking----------------------------------
c lines (i=constant,j=constant)
c 
c send to 19
c
      if (bd(19).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,ey,sz),1,linetype(1),neighbor(19),
     1                 6,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 14
c
      if (bd(14).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,sy-1,sz),1,linetype(1),neighbor(14),
     1                 6,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 23
c
      if (bd(23).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,sy,sz),1,linetype(1),neighbor(23),
     1                 7,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 10
c
      if (bd(10).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,ey+1,sz),1,linetype(1),neighbor(10),
     1                 7,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 14
c
      if (bd(14).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,sz),1,linetype(1),neighbor(14),
     1                 8,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 19
c
      if (bd(19).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,ey+1,sz),1,linetype(1),neighbor(19),
     1                 8,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 10
c
      if (bd(10).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,ey,sz),1,linetype(1),neighbor(10),
     1                 9,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 23
c
      if (bd(23).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,sy-1,sz),1,linetype(1),neighbor(23),
     1                 9,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c------------------------------------------------------------------------
c lines (i=constant,k=constant)
c 
c send to 25
c
      if (bd(25).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,sy,ez),1,linetype(2),neighbor(25),
     1                 10,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 12
c
      if (bd(12).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,sy,sz-1),1,linetype(2),neighbor(12),
     1                 10,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 21
c
      if (bd(21).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,sy,sz),1,linetype(2),neighbor(21),
     1                 11,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 16
c
      if (bd(16).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,sy,ez+1),1,linetype(2),neighbor(16),
     1                 11,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 12
c
      if (bd(12).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,sz),1,linetype(2),neighbor(12),
     1                 12,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 25
c
      if (bd(25).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,sy,ez+1),1,linetype(2),neighbor(25),
     1                 12,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 16
c
      if (bd(16).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,ez),1,linetype(2),neighbor(16),
     1                 13,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 21
c
      if (bd(21).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,sy,sz-1),1,linetype(2),neighbor(21),
     1                 13,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c------------------------------------------------------------------------
c lines (j=constant,k=constant)
c 
c send to 8
c
      if (bd(8).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,ey,ez),1,linetype(3),neighbor(8),
     1                 14,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 4
c
      if (bd(4).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx,sy-1,sz-1),1,linetype(3),neighbor(4),
     1                 14,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 2
c
      if (bd(2).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,ey,sz),1,linetype(3),neighbor(2),
     1                 15,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 6
c
      if (bd(6).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx,sy-1,ez+1),1,linetype(3),neighbor(6),
     1                 15,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 4
c
      if (bd(4).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,sz),1,linetype(3),neighbor(4),
     1                 16,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 8
c
      if (bd(8).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx,ey+1,ez+1),1,linetype(3),neighbor(8),
     1                 16,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 6
c
      if (bd(6).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,ez),1,linetype(3),neighbor(6),
     1                 17,comm3d,req(ireq),ierr)
# if cdebug
        tisend(2,1)=tisend(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 2
c
      if (bd(2).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx,ey+1,sz-1),1,linetype(3),neighbor(2),
     1                 17,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(2,1)=tirecv(2,1)+MPI_WTIME()-tmpi
# endif
      end if
c------------------------------------------------------------------------
# if cdebug
      nc=12-(bd(2)+bd(4)+bd(6)+bd(8)+bd(10)+bd(12)+bd(14)+bd(16)
     1      +bd(19)+bd(21)+bd(23)+bd(25))
      nisend(2,1)=nisend(2,1)+nc
      nirecv(2,1)=nirecv(2,1)+nc
# endif
# else
c----------------------------blocking------------------------------------
c lines (i=constant,j=constant)
c 
c send to 19 and receive from 14
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,ey,sz),1,linetype(1),neighbor(19),0,
     1                  a(sx-1,sy-1,sz),1,linetype(1),neighbor(14),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 23 and receive from 10
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,sy,sz),1,linetype(1),neighbor(23),1,
     1                  a(sx-1,ey+1,sz),1,linetype(1),neighbor(10),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 14 and receive from 19
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,sz),1,linetype(1),neighbor(14),0,
     1                  a(ex+1,ey+1,sz),1,linetype(1),neighbor(19),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 10 and receive from 23
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,ey,sz),1,linetype(1),neighbor(10),1,
     1                  a(ex+1,sy-1,sz),1,linetype(1),neighbor(23),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c------------------------------------------------------------------------
c lines (i=constant,k=constant)
c 
c send to 25 and receive from 12
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,sy,ez),1,linetype(2),neighbor(25),0,
     1                  a(sx-1,sy,sz-1),1,linetype(2),neighbor(12),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 21 and receive from 16
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,sy,sz),1,linetype(2),neighbor(21),1,
     1                  a(sx-1,sy,ez+1),1,linetype(2),neighbor(16),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 12 and receive from 25
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,sz),1,linetype(2),neighbor(12),0,
     1                  a(ex+1,sy,ez+1),1,linetype(2),neighbor(25),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 16 and receive from 21
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,ez),1,linetype(2),neighbor(16),1,
     1                  a(ex+1,sy,sz-1),1,linetype(2),neighbor(21),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c------------------------------------------------------------------------
c lines (j=constant,k=constant)
c 
c send to 8 and receive from 4
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,ey,ez),1,linetype(3),neighbor(8),0,
     1                  a(sx,sy-1,sz-1),1,linetype(3),neighbor(4),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 2 and receive from 6
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,ey,sz),1,linetype(3),neighbor(2),1,
     1                  a(sx,sy-1,ez+1),1,linetype(3),neighbor(6),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 4 and receive from 8
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,sz),1,linetype(3),neighbor(4),0,
     1                  a(sx,ey+1,ez+1),1,linetype(3),neighbor(8),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
c
c send to 6 and receive from 2
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,ez),1,linetype(3),neighbor(6),1,
     1                  a(sx,ey+1,sz-1),1,linetype(3),neighbor(2),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(2,1)=tsendrecv(2,1)+MPI_WTIME()-tmpi
# endif
# if cdebug
      nsendrecv(2,1)=nsendrecv(2,1)+12
# endif
# endif
c
# if cdebug
      timing(59)=timing(59)+MPI_WTIME()-tinitial
# endif
      return
      end
