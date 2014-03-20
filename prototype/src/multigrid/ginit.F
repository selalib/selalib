      subroutine ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi,pi,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,IOUT
      REALN p(sx-1:ex+1,sy-1:ey+1),r(sx-1:ex+1,sy-1:ey+1),
     1      f(sx-1:ex+1,sy-1:ey+1),hxi,hyi,wk,pi
c-----------------------------------------------------------------------
c Initialize the pressure, density, and right-hand side of the
c elliptic equation div(1/r*grad(p))=f
c
c Code      : tmgd2
c Called in : main
c Calls     : --
c-----------------------------------------------------------------------
      integer i,j
      REALN cnst,cx,cy,xi,yj
c
      do j=sy-1,ey+1
        do i=sx-1,ex+1
          p(i,j)=0.0d0
          r(i,j)=1.0d0
          f(i,j)=0.0d0
        end do
      end do
      cnst=-8.0d0*(pi*wk)**2
      cx=2.0d0*pi*wk
      cy=2.0d0*pi*wk
      do j=sy,ey
        yj=(float(j)-1.5d0)/hyi
        do i=sx,ex
          xi=(float(i)-1.5d0)/hxi
          f(i,j)=cnst*sin(cx*xi)*sin(cy*yj)
        end do
      end do
c
      return
      end
