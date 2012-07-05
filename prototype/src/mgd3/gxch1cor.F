      subroutine gxch1cor(sx,ex,sy,ey,sz,ez,a,comm3d,neighbor,
     1                    bd,cornertype,req,ireq,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,sz,ez,comm3d
      integer neighbor(26),bd(26),cornertype,req(52),ireq,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
c------------------------------------------------------------------------
c Subroutine to exchange 1*1*1 corners (i.e. points) of boundary data 
c between 'myid' and the 8 processes (i+,j+,k+), (i+,j+,k-), 
c (i+,j-,k+), (i+,j-,k-), (i-,j+,k+), (i-,j+,k-), (i-,j-,k+),
c (i-,j-,k-). Can be used to communicate scalar (p, r,...) as well as 
c vector quantities ((u,v,w), (ut,vt,wt)). The appropriate datatypes 
c have been defined in the 'type_mpi' subroutine.
c
c
c plane i
c       |           | 
c      ---------------                k
c       |           |               ^ 
c       |           |               |
c       |    myid   |               |
c       |           |              ------>
c       |           |               |     j
c      ---------------
c       |           | 
c
c plane i-1                  plane i+1
c     15|           |17          24|           |26 
c      ---------------            --------------- 
c       |           |              |           |
c       |           |              |           |
c       |           |              |           |
c       |           |              |           |
c       |           |              |           |
c      ---------------            ---------------
c     13|           |11          22|           |20
c
c Code      : tmgd3, test program for 3D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : gbackin, setdens, gdrbsor, gprbsor, gveloc, mgdrestr, 
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
c
c send to 26
c
      if (bd(26).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,ey,ez),1,cornertype,neighbor(26),
     1                 18,comm3d,req(ireq),ierr)
# if cdebug
        tisend(3,1)=tisend(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 13
c
      if (bd(13).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,sy-1,sz-1),1,cornertype,neighbor(13),
     1                 18,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(3,1)=tirecv(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 20
c
      if (bd(20).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,ey,sz),1,cornertype,neighbor(20),
     1                 19,comm3d,req(ireq),ierr)
# if cdebug
        tisend(3,1)=tisend(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 15
c
      if (bd(15).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,sy-1,ez+1),1,cornertype,neighbor(15),
     1                 19,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(3,1)=tirecv(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 24
c
      if (bd(24).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,sy,ez),1,cornertype,neighbor(24),
     1                 20,comm3d,req(ireq),ierr)
# if cdebug
        tisend(3,1)=tisend(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 11
c
      if (bd(11).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,ey+1,sz-1),1,cornertype,neighbor(11),
     1                 20,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(3,1)=tirecv(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 22
c
      if (bd(22).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,sy,sz),1,cornertype,neighbor(22),
     1                 21,comm3d,req(ireq),ierr)
# if cdebug
        tisend(3,1)=tisend(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 17
c
      if (bd(17).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,ey+1,ez+1),1,cornertype,neighbor(17),
     1                 21,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(3,1)=tirecv(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 13
c
      if (bd(13).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,sz),1,cornertype,neighbor(13),
     1                 22,comm3d,req(ireq),ierr)
# if cdebug
        tisend(3,1)=tisend(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 26
c
      if (bd(26).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,ey+1,ez+1),1,cornertype,neighbor(26),
     1                 22,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(3,1)=tirecv(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 15
c
      if (bd(15).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,ez),1,cornertype,neighbor(15),
     1                 23,comm3d,req(ireq),ierr)
# if cdebug
        tisend(3,1)=tisend(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 20
c
      if (bd(20).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,ey+1,sz-1),1,cornertype,neighbor(20),
     1                 23,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(3,1)=tirecv(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 11
c
      if (bd(11).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,ey,sz),1,cornertype,neighbor(11),
     1                 24,comm3d,req(ireq),ierr)
# if cdebug
        tisend(3,1)=tisend(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 24
c
      if (bd(24).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,sy-1,ez+1),1,cornertype,neighbor(24),
     1                 24,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(3,1)=tirecv(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 17
c
      if (bd(17).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,ey,ez),1,cornertype,neighbor(17),
     1                 25,comm3d,req(ireq),ierr)
# if cdebug
        tisend(3,1)=tisend(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 22
c
      if (bd(22).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,sy-1,sz-1),1,cornertype,neighbor(22),
     1                 25,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(3,1)=tirecv(3,1)+MPI_WTIME()-tmpi
# endif
      end if
c------------------------------------------------------------------------
# if cdebug
      nc=8-(bd(26)+bd(13)+bd(20)+bd(15)+bd(24)+bd(11)+bd(22)+bd(17))
      nisend(3,1)=nisend(3,1)+nc
      nirecv(3,1)=nirecv(3,1)+nc
# endif
# else
c----------------------------blocking------------------------------------
c send to 26 and receive from 13
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,ey,ez),1,cornertype,neighbor(26),0,
     1                  a(sx-1,sy-1,sz-1),1,cornertype,neighbor(13),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(3,1)=tsendrecv(3,1)+MPI_WTIME()-tmpi
# endif
c
c send to 20 and receive from 15
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,ey,sz),1,cornertype,neighbor(20),1,
     1                  a(sx-1,sy-1,ez+1),1,cornertype,neighbor(15),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(3,1)=tsendrecv(3,1)+MPI_WTIME()-tmpi
# endif
c
c send to 24 and receive from 11
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,sy,ez),1,cornertype,neighbor(24),0,
     1                  a(sx-1,ey+1,sz-1),1,cornertype,neighbor(11),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(3,1)=tsendrecv(3,1)+MPI_WTIME()-tmpi
# endif
c
c send to 22 and receive from 17
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,sy,sz),1,cornertype,neighbor(22),1,
     1                  a(sx-1,ey+1,ez+1),1,cornertype,neighbor(17),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(3,1)=tsendrecv(3,1)+MPI_WTIME()-tmpi
# endif
c
c send to 13
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,sz),1,cornertype,neighbor(13),0,
     1                  a(ex+1,ey+1,ez+1),1,cornertype,neighbor(26),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(3,1)=tsendrecv(3,1)+MPI_WTIME()-tmpi
# endif
c
c send to 15 and receive from 20
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,ez),1,cornertype,neighbor(15),1,
     1                  a(ex+1,ey+1,sz-1),1,cornertype,neighbor(20),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(3,1)=tsendrecv(3,1)+MPI_WTIME()-tmpi
# endif
c
c send to 11 and receive from 24
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,ey,sz),1,cornertype,neighbor(11),0,
     1                  a(ex+1,sy-1,ez+1),1,cornertype,neighbor(24),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(3,1)=tsendrecv(3,1)+MPI_WTIME()-tmpi
# endif
c
c send to 17 and receive from 22
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,ey,ez),1,cornertype,neighbor(17),1,
     1                  a(ex+1,sy-1,sz-1),1,cornertype,neighbor(22),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(3,1)=tsendrecv(3,1)+MPI_WTIME()-tmpi
# endif
# if cdebug
      nsendrecv(3,1)=nsendrecv(3,1)+8
# endif
# endif
c
# if cdebug
      timing(60)=timing(60)+MPI_WTIME()-tinitial
# endif
      return
      end
