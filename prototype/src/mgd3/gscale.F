      subroutine gscale(sx,ex,sy,ey,sz,ez,a,avo,acorr,comm3d,nx,ny,nz,
     1                  isol,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,sz,ez,nx,ny,nz,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),avo,acorr
      integer comm3d,isol
c------------------------------------------------------------------------
c Rescale the field a so that its average inside the domain
c remains constant and equal to avo. For the density,avo should
c be rro, this ensures conservation of mass. For the pressure,
c avo should be 0 so that the average pressure does not drift
c away from 0, which is the initial value. For the density, also
c eliminate the numerical overshoots and undershoots in the solution,
c so that can consider higher density ratios; this feature is active
c only if the compiler directive UNDERSHOOT is set to 1 in the
c "compdir.inc" file.
c
c Code      : tmgd3, test program for 3D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver, gdrbsor, gprbsor
c Calls     : MPI_ALLREDUCE
c------------------------------------------------------------------------
      REALN avloc,av
      integer i,j,k,ierr
# if cdebug
      double precision tinitial,tmpi
      tinitial=MPI_WTIME()
# endif
c
c determine average value
c
      avloc=0.0d0
      do k=sz,ez
        do j=sy,ey
          do i=sx,ex
            avloc=avloc+a(i,j,k)
          end do
        end do
      end do
c
c global reduce across all process
c
# if cdebug
      tmpi=MPI_WTIME()
# endif
# if double_precision
      call MPI_ALLREDUCE(avloc,av,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     1                   comm3d,ierr)
# else
      call MPI_ALLREDUCE(avloc,av,1,MPI_REAL,MPI_SUM,comm3d,ierr)
# endif
# if cdebug
      nallreduce=nallreduce+1
      tallreduce=tallreduce+MPI_WTIME()-tmpi
# endif
      av=av/float(nx*ny*nz)
c
c do correction
c
      acorr=avo-av
      do k=sz,ez
        do j=sy,ey
          do i=sx,ex
            a(i,j,k)=a(i,j,k)+acorr
          end do
        end do
      end do
c
# if cdebug
      timing(49)=timing(49)+MPI_WTIME()-tinitial
# endif
      return
      end 

