      subroutine gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,pi,nx,ny,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,comm2d,IOUT,nx,ny
      REALN p(sx-1:ex+1,sy-1:ey+1),wk,hxi,hyi,pi
c-----------------------------------------------------------------------
c Calculate the error between the numerical and exact solution to
c the test problem.
c
c Code      : tmgd2
c Called in : main
c Calls     : MPI_ALLREDUCE
c-----------------------------------------------------------------------
      integer i,j,ierr
      REALN errloc,err,cx,cy,exact
c
c calculate local error
c
      cx=2.0d0*pi*wk
      cy=2.0d0*pi*wk
      errloc=0.0d0
      do j=sy,ey
        yj=(float(j)-1.5d0)/hyi
        do i=sx,ex
          xi=(float(i)-1.5d0)/hxi
          exact=sin(cx*xi)*sin(cy*yj)
          errloc=errloc+abs(p(i,j)-exact)
        end do
      end do
c
c calculate global error
c
# if double_precision
      call MPI_ALLREDUCE(errloc,err,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     1                   comm2d,ierr)
# else
      call MPI_ALLREDUCE(errloc,err,1,MPI_REAL,MPI_SUM,comm2d,ierr)
# endif
      write(IOUT,100) errloc/float(nx*ny),err/float(nx*ny)
100   format(/,'Local error: ',e13.6,'  total error: ',e13.6,/)
c
      return
      end
