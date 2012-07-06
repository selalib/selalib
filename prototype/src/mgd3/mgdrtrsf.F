      subroutine mgdrtrsf(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,rc,
     1                    sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,rf,
     2                    comm3dp,myid,neighbor,bd,planetype,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc
      integer sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf
      integer comm3dp,myid,neighbor(26),bd(26),planetype(3),IOUT
      REALN rc(sxc-1:exc+1,syc-1:eyc+1,szc-1:ezc+1)
      REALN rf(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
c------------------------------------------------------------------------
c Transfer values of the density from a finer to a coarser grid level.
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : gxch1pla, MPI_WAITALL (non-blocking version)
c------------------------------------------------------------------------
      integer i,j,k,ic,jc,kc,iinc,jinc,kinc,i1,i2,j1,j2,k1,k2
      integer ireq,req(52)
# if NBLOCKGR
      integer status(MPI_STATUS_SIZE,52),ierr
# endif
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      if (nxc.lt.nxf) then
        iinc=2
        i1=1
        i2=0
      else 
        iinc=1 
        i1=0
        i2=1
      end if
      if (nyc.lt.nyf) then
        jinc=2
        j1=1
        j2=0
      else
        jinc=1
        j1=0
        j2=1
      end if
      if (nzc.lt.nzf) then
        kinc=2
        k1=1
        k2=0
      else
        kinc=1
        k1=0
        k2=1
      end if
      do kc=szc,ezc
        k=k1*(2*kc-1)+k2*kc
        do jc=syc,eyc
          j=j1*(2*jc-1)+j2*jc
          do ic=sxc,exc
            i=i1*(2*ic-1)+i2*ic
            rc(ic,jc,kc)=rf(i,j,k)
          end do
        end do
      end do
c
c exchange the plane boundary values at coarse level
c
# if NBLOCKGR
      ireq=0
# endif
      call gxch1pla(sxc,exc,syc,eyc,szc,ezc,rc,comm3dp,neighbor,
     1              bd,planetype,req,ireq,IOUT)
# if NBLOCKGR
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
      nwaitall=nwaitall+1
      twaitall=twaitall+MPI_WTIME()-tmpi
# endif
# endif
c
# if cdebug
      timing(86)=timing(86)+MPI_WTIME()-tinitial
# endif
      return
      end
