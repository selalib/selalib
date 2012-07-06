      subroutine ginit(sx,ex,sy,ey,sz,ez,p,r,f,wk,hxi,hyi,hzi,pi,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,sz,ez,IOUT
      REALN p(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),
     1      r(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),
     2      f(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),hxi,hyi,hzi,wk,pi
c-----------------------------------------------------------------------
c Initialize the pressure, density, and right-hand side of the
c elliptic equation div(1/r*grad(p))=f
c
c Code      : tmgd3, test program for 3D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : --
c-----------------------------------------------------------------------
      integer i,j,k
      REALN cnst,cx,cy,cz,xi,yj,zk
c
      do k=sz-1,ez+1
        do j=sy-1,ey+1
          do i=sx-1,ex+1
            p(i,j,k)=0.0d0
            r(i,j,k)=1.0d0
            f(i,j,k)=0.0d0
          end do
        end do
      end do
      cnst=-12.0d0*(pi*wk)**2
      cx=2.0d0*pi*wk
      cy=2.0d0*pi*wk
      cz=2.0d0*pi*wk
      do k=sz,ez
        zk=(float(k)-1.5d0)/hzi
        do j=sy,ey
          yj=(float(j)-1.5d0)/hyi
          do i=sx,ex
            xi=(float(i)-1.5d0)/hxi
            f(i,j,k)=cnst*sin(cx*xi)*sin(cy*yj)*sin(cz*zk)
          end do
        end do
      end do
c
      return
      end
