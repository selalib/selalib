      subroutine gerr(sx,ex,sy,ey,sz,ez,p,comm3d,wk,hxi,hyi,hzi,
     1                pi,nx,ny,nz,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,sz,ez,comm3d,IOUT,nx,ny,nz
      REALN p(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),wk,hxi,hyi,hzi,pi
c-----------------------------------------------------------------------
c Calculate the error between the numerical and exact solution to
c the test problem.
c
c Code      : tmgd3, test program for 3D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : MPI_ALLREDUCE
c-----------------------------------------------------------------------
      integer i,j,k,ierr
      REALN errloc,err,cx,cy,cz,exact,zk,yj,xi
c
c calculate local error
c
      cx=2.0d0*pi*wk
      cy=2.0d0*pi*wk
      cz=2.0d0*pi*wk
      errloc=0.0d0
      do k=sz,ez
        zk=(float(k)-1.5d0)/hzi
        do j=sy,ey
          yj=(float(j)-1.5d0)/hyi
          do i=sx,ex
            xi=(float(i)-1.5d0)/hxi
            exact=sin(cx*xi)*sin(cy*yj)*sin(cz*zk)
            errloc=errloc+abs(p(i,j,k)-exact)
          end do
        end do
      end do
c
c calculate global error
c
      call MPI_ALLREDUCE(errloc,err,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     1                   comm3d,ierr)
      write(IOUT,100) errloc/float(nx*ny*nz),err/float(nx*ny*nz)
100   format(/,'Local error: ',e13.6,'  total error: ',e13.6,/)
c
      return
      end
