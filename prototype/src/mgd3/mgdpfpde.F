      subroutine mgdpfpde(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,cof,r,
     1                    xl,yl,zl,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,IOUT
      REALN cof(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1,8)
      REALN r(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1),xl,yl,zl
c------------------------------------------------------------------------
c Determine coefficients for the pressure equation at the finest grid
c level. These coefficients involve densities half-way between the
c pressure and density nodes.
c
c cof array: 1 -> i-1
c            2 -> i+1
c            3 -> j-1
c            4 -> j+1
c            5 -> k-1
c            6 -> k+1
c            7 -> central
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      REALN dlx,todlxx,dly,todlyy,dlz,todlzz,rijk,c1,c2,c3,c4,c5,c6
      integer i,j,k
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      dlx=xl/float(nxf-1)
      todlxx=2.0d0/(dlx*dlx)
      dly=yl/float(nyf-1)
      todlyy=2.0d0/(dly*dly)
      dlz=zl/float(nzf-1)
      todlzz=2.0d0/(dlz*dlz)
      do k=szf,ezf
        do j=syf,eyf
          do i=sxf,exf
            rijk=r(i,j,k)
            c1=todlxx/(rijk+r(i-1,j,k))
            c2=todlxx/(rijk+r(i+1,j,k))
            c3=todlyy/(rijk+r(i,j-1,k))
            c4=todlyy/(rijk+r(i,j+1,k))
            c5=todlzz/(rijk+r(i,j,k-1))
            c6=todlzz/(rijk+r(i,j,k+1))
            cof(i,j,k,1)=c1
            cof(i,j,k,2)=c2
            cof(i,j,k,3)=c3
            cof(i,j,k,4)=c4
            cof(i,j,k,5)=c5
            cof(i,j,k,6)=c6
            cof(i,j,k,7)=-(c1+c2+c3+c4+c5+c6)
          end do
        end do
      end do
c
# if cdebug
      timing(84)=timing(84)+MPI_WTIME()-tinitial
# endif
      return
      end
