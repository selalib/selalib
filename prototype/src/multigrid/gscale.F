      subroutine gscale(sx,ex,sy,ey,a,avo,acorr,comm2d,nx,ny,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,nx,ny,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1),avo,acorr
      integer comm2d
c------------------------------------------------------------------------
c Rescale the field a so that its average inside the domain
c remains constant and equal to avo. For the density,avo should
c be rro, this ensures conservation of mass. For the pressure,
c avo should be 0 so that the average pressure does not drift
c away from 0, which is the initial value.
c
c Code      : tmgd2
c Called in : mgdsolver
c Calls     : MPI_ALLREDUCE
c------------------------------------------------------------------------
      REALN avloc,av
      integer i,j,ierr
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c determine average value
c
      avloc=0.0d0
      do j=sy,ey
        do i=sx,ex
          avloc=avloc+a(i,j)
        end do
      end do
c
c global reduce across all process
c
# if double_precision
      call MPI_ALLREDUCE(avloc,av,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     1                   comm2d,ierr)
# else
      call MPI_ALLREDUCE(avloc,av,1,MPI_REAL,MPI_SUM,comm2d,ierr)
# endif
# if cdebug
      nallreduce=nallreduce+1
# endif
      av=av/float(nx*ny)
c
c do correction
c
      acorr=avo-av
      do j=sy,ey
        do i=sx,ex
          a(i,j)=a(i,j)+acorr
        end do
      end do
c
# if cdebug
      timing(49)=timing(49)+MPI_WTIME()-tinitial
# endif
      return
      end 

