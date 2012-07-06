      subroutine gxch1pla(sx,ex,sy,ey,sz,ez,a,comm3d,neighbor,
     1                    bd,planetype,req,ireq,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,sz,ez,comm3d
      integer neighbor(26),bd(26),planetype(3),req(52),ireq,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
c------------------------------------------------------------------------
c Subroutine to exchange planes of thickness 1 of boundary data between
c 'myid' and the 6 processes in the directions i+, i-, j+, j-, k+, k-. 
c Can be used to communicate scalar (p, r,...) as well as vector 
c quantities ((u,v,w), (ut,vt,wt)). The appropriate datatypes have 
c been defined in the 'type_mpi' subroutine.
c
c plane i
c       |     7     |
c      ---------------                k
c       |           |               ^ 
c       |           |               |
c      5|    myid   |1              |
c       |           |              ------>
c       |           |               |     j
c      ---------------
c       |     3     |
c
c plane i-1                  plane i+1
c       |           |              |           | 
c      ---------------            --------------- 
c       |           |              |           |
c       |           |              |           |
c       |     9     |              |    18     |
c       |           |              |           |
c       |           |              |           |
c      ---------------            ---------------
c       |           |              |           |
c
c Code      : tmgd3, test program for 3D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : setdens,gbackin, gadvect, gdrbsor, gprbsor, gveloc, 
c             mgdrelax, mgdrestr, mgdrtrsf, 
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
c planes i=constant
c
c send to 18
c
      if (bd(18).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(ex,sy,sz),1,planetype(1),neighbor(18),
     1                 0,comm3d,req(ireq),ierr)
# if cdebug
        tisend(1,1)=tisend(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 9
c
      if (bd(9).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx-1,sy,sz),1,planetype(1),neighbor(9),
     1                 0,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(1,1)=tirecv(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 9
c
      if (bd(9).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,sz),1,planetype(1),neighbor(9),
     1                 1,comm3d,req(ireq),ierr)
# if cdebug
        tisend(1,1)=tisend(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 18
c
      if (bd(18).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(ex+1,sy,sz),1,planetype(1),neighbor(18),
     1                 1,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(1,1)=tirecv(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c------------------------------------------------------------------------
c planes j=constant
c
c send to 1
c
      if (bd(1).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,ey,sz),1,planetype(2),neighbor(1),
     1                 2,comm3d,req(ireq),ierr)
# if cdebug
        tisend(1,1)=tisend(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 5
c
      if (bd(5).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx,sy-1,sz),1,planetype(2),neighbor(5),
     1                 2,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(1,1)=tirecv(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 5
c
      if (bd(5).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,sz),1,planetype(2),neighbor(5),
     1                 3,comm3d,req(ireq),ierr)
# if cdebug
        tisend(1,1)=tisend(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 1
c
      if (bd(1).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx,ey+1,sz),1,planetype(2),neighbor(1),
     1                 3,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(1,1)=tirecv(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c------------------------------------------------------------------------
c planes k=constant
c
c send to 7
c
      if (bd(7).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,ez),1,planetype(3),neighbor(7),
     1                 4,comm3d,req(ireq),ierr)
# if cdebug
        tisend(1,1)=tisend(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 3
c
      if (bd(3).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx,sy,sz-1),1,planetype(3),neighbor(3),
     1                 4,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(1,1)=tirecv(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c send to 3
c
      if (bd(3).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_ISEND(a(sx,sy,sz),1,planetype(3),neighbor(3),
     1                 5,comm3d,req(ireq),ierr)
# if cdebug
        tisend(1,1)=tisend(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c
c receive from 7
c
      if (bd(7).eq.0) then
        ireq=ireq+1
# if cdebug
        tmpi=MPI_WTIME()
# endif
        call MPI_IRECV(a(sx,sy,ez+1),1,planetype(3),neighbor(7),
     1                 5,comm3d,req(ireq),ierr)
# if cdebug
        tirecv(1,1)=tirecv(1,1)+MPI_WTIME()-tmpi
# endif
      end if
c------------------------------------------------------------------------
# if cdebug
      nc=6-(bd(1)+bd(3)+bd(5)+bd(7)+bd(9)+bd(18))
      nisend(1,1)=nisend(1,1)+nc
      nirecv(1,1)=nirecv(1,1)+nc
# endif
# else
c----------------------------blocking------------------------------------
c planes i=constant
c
c send to 18 and receive from 9
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(ex,sy,sz),1,planetype(1),neighbor(18),0,
     1                  a(sx-1,sy,sz),1,planetype(1),neighbor(9),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(1,1)=tsendrecv(1,1)+MPI_WTIME()-tmpi
# endif
c
c send to 9 and receive from 18
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,sz),1,planetype(1),neighbor(9),1,
     1                  a(ex+1,sy,sz),1,planetype(1),neighbor(18),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(1,1)=tsendrecv(1,1)+MPI_WTIME()-tmpi
# endif
c------------------------------------------------------------------------
c planes j=constant
c
c send to 1 and receive from 5
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,ey,sz),1,planetype(2),neighbor(1),0,
     1                  a(sx,sy-1,sz),1,planetype(2),neighbor(5),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(1,1)=tsendrecv(1,1)+MPI_WTIME()-tmpi
# endif
c
c send to 5 and receive from 1
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,sz),1,planetype(2),neighbor(5),1,
     1                  a(sx,ey+1,sz),1,planetype(2),neighbor(1),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(1,1)=tsendrecv(1,1)+MPI_WTIME()-tmpi
# endif
c------------------------------------------------------------------------
c planes k=constant
c
c send to 7 and receive from 3
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,ez),1,planetype(3),neighbor(7),0,
     1                  a(sx,sy,sz-1),1,planetype(3),neighbor(3),0,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(1,1)=tsendrecv(1,1)+MPI_WTIME()-tmpi
# endif
c
c send to 3 and receive from 7
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_SENDRECV(a(sx,sy,sz),1,planetype(3),neighbor(3),1,
     1                  a(sx,sy,ez+1),1,planetype(3),neighbor(7),1,
     2                  comm3d,status,ierr)
# if cdebug
      tsendrecv(1,1)=tsendrecv(1,1)+MPI_WTIME()-tmpi
# endif
# if cdebug
      nsendrecv(1,1)=nsendrecv(1,1)+6
# endif
# endif
c
# if cdebug
      timing(58)=timing(58)+MPI_WTIME()-tinitial
# endif
      return
      end
